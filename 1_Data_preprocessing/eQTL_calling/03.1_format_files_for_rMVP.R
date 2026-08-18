## Format PANZEA and WIDIV combined markers for mapping (GWAS and eQTL) with rMVP
library(data.table)
library(tidyverse)
library(rMVP)

## Set paths to geno files
my_geno_path = path.expand("~/eqtl_pdi_coexpr/data/02_genotypic_data/01_individual_beds/")

## List da bed files
setwd(my_geno_path)
my_geno = list.files(pattern = "\\.bed$")

## Set path to expression files
my_phenos_path = path.expand("~/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/")
my_phenos = list.files(my_phenos_path, pattern = "\\.expression.bed.gz$")

setwd("~/eqtl_pdi_coexpr/data/07_rMVP_data/")
for (i in seq_along(my_phenos)) {
  mf_name = gsub(".expression.bed.gz","",my_phenos[i])
  
  ge = fread(paste0(my_phenos_path,my_phenos[i])) %>% 
    select(-c("#chr",start,end))
  gene_ids = ge$gene_id
  ge = t(ge[,-1])
  taxa_ids = rownames(ge)
  ge = as_tibble(ge) %>% 
    setNames(., c(gene_ids)) %>% 
    mutate(Taxa = taxa_ids, .before = 1)
  
  ## Write the file
  fwrite(ge, paste0(mf_name,".txt") , sep = "\t", col.names = TRUE, row.names = FALSE)
  my_phenos_tmp = paste0(mf_name,".txt")
  
  ## Make the datasets 
  MVP.Data(fileBed = paste0(my_geno_path,mf_name),
           fileKin=TRUE,
           filePhe=my_phenos_tmp,
           filePC=TRUE,
           SNP.impute = "Major",
           ncpus = 20,
           out=mf_name)

}

print("DONE")
# system("rm *.txt")

