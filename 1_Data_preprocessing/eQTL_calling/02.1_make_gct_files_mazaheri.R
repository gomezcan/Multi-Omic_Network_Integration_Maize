## Format expression datasets from Peng to be compatible for TensorQTL, and FarmCPU
library(data.table)
library(tidyverse)

## For TensorQTL need to have convert data to .gct file 
setwd("~/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/")

## Make .gct file for WIDIV 942 genos but keep only 602 lines with phenolics
widiv_942_rTPM = fread("widiv942_rTPM.tsv")
widiv_942_Reads = fread("widiv942_ReadCount.tsv")
widiv_942_meta =  fread("01.meta_widiv942.txt", na.strings = "")

## Isolate the different sets of genotypes

## 602 lines with phenolic infomation
widiv_602_meta = as_tibble(widiv_942_meta) %>% 
  select(SampleID, Genotype, phenolics) %>% 
  filter(complete.cases(.))

widiv_602_lines = as_tibble(widiv_942_meta) %>% 
  select(SampleID, Genotype, phenolics) %>% 
  filter(complete.cases(.)) %>% 
  pull(Genotype)

## 304 lines with phenolic infomation and WGS SNP data
widiv_304_meta = as_tibble(widiv_942_meta) %>% 
  select(SampleID, Genotype, phenolics, hapmap) %>% 
  filter(complete.cases(.))

widiv_304_lines = as_tibble(widiv_942_meta) %>% 
  select(SampleID, phenolics, hapmap) %>% 
  filter(complete.cases(.)) %>% 
  pull(hapmap)

## 187 lines with phenolic infomation, WGS SNP data, and expression data in more than one tissue
setwd("~/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/")
lines_krem = fread("0x_geno_names_across_studies.txt", na.strings = "") %>% filter(complete.cases(.))

widiv_187_meta = as_tibble(widiv_942_meta) %>% 
  select(SampleID, Genotype, phenolics, hapmap) %>% 
  filter(hapmap%in%lines_krem$genotype_hapmap)

widiv_187_lines = as_tibble(widiv_942_meta) %>% 
  select(SampleID, Genotype, phenolics, hapmap) %>% 
  filter(hapmap%in%lines_krem$genotype_hapmap) %>% 
  pull(hapmap)

## Make function to create .gct file given metadata, number of genos included, expression matrix
################################################################################################
# my_expression = expression matrix
# my_meta = the meta information corresponding to the expression matrix
# my_lines = the lines to subset from the expression matrix and write a .gct file from
################################################################################################
my_gct_function = function(my_expression, my_meta, my_lines){
  e = my_expression %>% 
    select(gid, my_meta$SampleID)
  
  empt_row = rep(NA, length(my_lines))
  gene_names = e$gid
  
  e1 = as_tibble(e[,-1]) %>% 
    rbind(my_lines,.) %>% 
    cbind(as_tibble(c("Description",gene_names)),.) %>% 
    cbind(as_tibble(c("NAME",gene_names)),.) %>% 
    rbind(empt_row,.) %>% 
    rbind(empt_row,.)
  
  e1[1,1] = "#1.2"
  e1[2,1] = as.character(ncol(e)-1)
  e1[2,2] = as.character(nrow(e))
  
  t_names_write = as_tibble(my_lines) %>% 
    cbind(my_lines) %>% 
    setNames(., c("sample_id","participant_id"))
  
  return(list(e1, t_names_write))
}

x_187_raw = my_gct_function(my_expression = widiv_942_Reads, my_meta = widiv_187_meta, my_lines = widiv_187_lines)
x_304_raw = my_gct_function(my_expression = widiv_942_Reads, my_meta = widiv_304_meta, my_lines = widiv_304_lines)
x_602_raw = my_gct_function(my_expression = widiv_942_Reads, my_meta = widiv_602_meta, my_lines = widiv_602_lines)

x_187_tpm = my_gct_function(my_expression = widiv_942_rTPM, my_meta = widiv_187_meta, my_lines = widiv_187_lines)
x_304_tpm = my_gct_function(my_expression = widiv_942_rTPM, my_meta = widiv_304_meta, my_lines = widiv_304_lines)
x_602_tpm = my_gct_function(my_expression = widiv_942_rTPM, my_meta = widiv_602_meta, my_lines = widiv_602_lines)

setwd("~/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/")
fwrite(x_187_raw[[1]], "001.ReadCount_WIDIV_seedling_V1_expression_187g_agpv4.gct",quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(x_187_tpm[[1]], "001.rTPM_WIDIV_seedling_V1_expression_187g_agpv4.gct",quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(x_187_tpm[[2]], "001.taxa_lookup_187g_agpv4.txt",quote = FALSE, sep = "\t", col.names = TRUE)

fwrite(x_304_raw[[1]], "002.ReadCount_WIDIV_seedling_V1_expression_304g_agpv4.gct",quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(x_304_tpm[[1]], "002.rTPM_WIDIV_seedling_V1_expression_304g_agpv4.gct",quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(x_304_tpm[[2]], "002.taxa_lookup_304g_agpv4.txt",quote = FALSE, sep = "\t", col.names = TRUE)

fwrite(x_602_raw[[1]], "003.ReadCount_WIDIV_seedling_V1_expression_602g_agpv4.gct",quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(x_602_tpm[[1]], "003.rTPM_WIDIV_seedling_V1_expression_602g_agpv4.gct",quote = FALSE, sep = "\t", col.names = FALSE)
fwrite(x_602_tpm[[2]], "003.taxa_lookup_602g_agpv4.txt",quote = FALSE, sep = "\t", col.names = TRUE)


tmp = data.table::fread("001.taxa_lookup_187g_agpv4.txt")
