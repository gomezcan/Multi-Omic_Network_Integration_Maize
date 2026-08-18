## Format raw TPM expression datasets from Peng to be compatible for TensorQTL, and FarmCPU
## Will do TMM normalization following gtex pipeline

library(data.table)
library(tidyverse)
library(dplyr)
library(trqwe)

## Read in Kremling 2018 tissue*genotype dataset and meta data corresponding to names in combined SNP file prepared in previous step
setwd("~/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/")
x = mcreadRDS("cpm.rds")
meta_tmp = fread("meta.csv") %>% select(Genotype, Genotype_mazaheri,Genotype_hapmap_and_widiv) %>% distinct()
meta = x$th_m %>% left_join(., meta_tmp, by = c("Genotype"))

## Isolate the read counts then normalized expression matrix and keep the lines overlapping with widiv, hapmap, and phenolics
lines_keep = fread("0x_geno_names_across_studies.txt", na.strings = "") %>% filter(complete.cases(.))

## Do calculations to get normalized TPM across average of genotype replicates per tissue
my_x = x$tm_m %>% 
  left_join(., meta, by = "SampleID") %>% 
  filter(Genotype_hapmap_and_widiv%in%lines_keep$genotype_hapmap) %>%
  group_by(SampleID) %>% 
  mutate(gene_length = rCPM/rFPKM) %>% 
  mutate(RPK = ReadCount/gene_length) %>%
  mutate(scale_fact = sum(RPK,na.rm = TRUE)/1e6) %>% 
  mutate(rTPM = RPK/scale_fact)

## Work through one tissue at a time to write out data
setwd("~/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/")
my_tissues = unique(my_x$Tissue)

for(j in seq_along(my_tissues)){
  i = my_tissues[j]
  my_x1 = as_tibble(my_x) %>% 
    filter(Tissue==i) %>%
    pivot_wider(id_cols = c(gid), names_from = Genotype_hapmap_and_widiv, values_from = rTPM) %>% 
    replace(is.na(.), 0)
  
  my_x2 = as_tibble(my_x) %>% 
    filter(Tissue==i) %>% 
    pivot_wider(id_cols = c(gid), names_from = Genotype_hapmap_and_widiv, values_from = ReadCount) %>% 
    replace(is.na(.), 0)
  
  t_names1 = colnames(my_x1[2:ncol(my_x1)])
  g_names1 = my_x1$gid
  e_row1 = rep(NA, length(t_names1))
  
  my_df1 = as_tibble(my_x1[,-1]) %>%
    rbind(t_names1,.) %>% 
    cbind(as_tibble(c("Description",g_names1)),.) %>% 
    cbind(as_tibble(c("NAME",g_names1)),.) %>% 
    rbind(e_row1,.) %>% 
    rbind(e_row1,.)
  
  my_df1[1,1] = "#1.2"
  my_df1[2,1] = as.character(ncol(my_x1)-1)
  my_df1[2,2] = as.character(nrow(my_x1))
  
  t_names2 = colnames(my_x2[2:ncol(my_x2)])
  g_names2 = my_x2$gid
  e_row2 = rep(NA, length(t_names2))
  
  my_df2 = as_tibble(my_x2[,-1]) %>%
    rbind(t_names2,.) %>% 
    cbind(as_tibble(c("Description",g_names2)),.) %>% 
    cbind(as_tibble(c("NAME",g_names2)),.) %>% 
    rbind(e_row2,.) %>% 
    rbind(e_row2,.)
  my_df2[1,1] = "#1.2"
  my_df2[2,1] = as.character(ncol(my_x2)-1)
  my_df2[2,2] = as.character(nrow(my_x2))
  
  t_names_write = as_tibble(t_names1) %>% 
    cbind(t_names1) %>% 
    setNames(., c("sample_id","participant_id"))
  
  fwrite(my_df1, paste0("00",j,".","rTPM_Kremling_",i,"_expression_",length(t_names1),"g_agpv4.gct"),
         quote = FALSE, sep = "\t", col.names = FALSE)
  
  fwrite(my_df2, paste0("00",j,".","ReadCount_Kremling_",i,"_expression_",length(t_names2),"g_agpv4.gct"),
         quote = FALSE, sep = "\t", col.names = FALSE)
    
  fwrite(t_names_write, paste0("00",j,".","taxa_lookup_",i,"_",length(t_names2),"g_agpv4.txt"),
         quote = FALSE, sep = "\t", col.names = TRUE)
  
  print(paste0("Done with file  ",j, "   out of  ",length(my_tissues)))
  
}


