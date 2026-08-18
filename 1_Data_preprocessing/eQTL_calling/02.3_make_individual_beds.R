## Make bed files for each data set to be analyzed
library(data.table)
library(tidyverse)

setwd("~/eqtl_pdi_coexpr/data/03_expression_data/01_kremling_2018/01_gct_files/")
t_files_krem = list.files(pattern = "\\.txt$")

f_krem = list()
for(i in seq_along(t_files_krem)){
  f_krem[[i]] = fread(t_files_krem[i]) %>% 
    select(sample_id)
}

setwd("~/eqtl_pdi_coexpr/data/03_expression_data/02_mazaheri_2019/01_gct_files/")
t_files_maza = list.files(pattern = "\\.txt$")

f_maza = list()
for(i in seq_along(t_files_maza)){
  f_maza[[i]] = fread(t_files_maza[i]) %>% 
    select(sample_id)
}


setwd("~/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/")
for(i in seq_along(t_files_krem)){
  fwrite(f_krem[[i]], file = paste0(t_files_krem[i]), col.names = FALSE, row.names = FALSE, sep = "\t") 
}
for(i in seq_along(t_files_maza)){
  fwrite(f_maza[[i]], file = paste0(t_files_maza[i]), col.names = FALSE, row.names = FALSE, sep = "\t") 
}

