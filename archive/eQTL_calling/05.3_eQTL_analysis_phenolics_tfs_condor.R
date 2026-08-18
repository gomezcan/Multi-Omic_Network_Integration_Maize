## Create a copy of the directory where all files live
# cp -r eqtl_pdi_coexpr/ eqtl_pdi_coexpr_condor
#
# ## Install R packages in the R build on condor
# install.packages("devtools")
# devtools::install_github("xiaolei-lab/rMVP")
# devtools::install_github("tidyverse/tidyverse")
# devtools::install_github("Rdatatable/data.table")
# install.packages("tictoc")
setwd("~/eqtl_pdi_coexpr_condor/")

args = commandArgs(trailingOnly = TRUE)

# k = as.numeric(args[1])
# sp = as.numeric(args[2])
# chrom = as.numeric(args[3])
k = 9
sp = 1
chrom = 1

## arguments from the command line:
# k = the file name suffix - should only be 1 ... 10
# sp = the chunk being processed. This will come from the text file

## Permutations run for each dataset separately using tensorQTL
library(rMVP)
library(data.table)
library(tidyverse)
library(tictoc)

## Set paths to data in rMVP format
my_phenos_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/data/07_rMVP_data/"))
my_tf_genes_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/data/05_pdi_data/"))
my_permutations_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/results/02_eQTL_permutations/"))
my_results_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/results/04_eQTL_results_phenolics_and_tfs/"))

## Create list of files to iterate through
setwd(my_phenos_path)
my_files = list.files(pattern = "\\.ind$") %>% gsub(".geno.ind","",.)

# ## Read in list of phenolic genes and TFs if running only subset
# setwd(my_tf_genes_path)
# my_phenolic_genes = fread("PhenolicGenes_fabio.csv", header = TRUE) %>% select(gene)
# my_tf_genes = fread("Zm_TF_CoR_Full_list.txt", header = FALSE) %>% setNames(., "gene")
# my_genes = bind_rows(my_tf_genes,my_phenolic_genes)

## Read in snp data
setwd(my_phenos_path)
myGD = attach.big.matrix(paste0(my_phenos_path,my_files,".geno.desc"))

## Read in map data
map = fread(paste0(my_phenos_path,my_files,".geno.map"))

## Keep daters for chromosome being analyzed
map_idx = which(map$CHROM == paste0("chr",all_of(chrom)))
myGD = deepcopy(myGD,rows = map_idx)
map = map %>% filter(CHROM == paste0("chr",all_of(chrom)))

## Read in Kinship
kinship = attach.big.matrix(paste0(my_phenos_path,my_files,".kin.desc"))

## Read in expression data and filter for phenolic and TF genes present
phenos = fread(paste0(my_phenos_path,my_files,".phe"))
# genes_in_expression_dat = my_genes %>% 
#   filter(gene%in%colnames(phenos)) %>% 
#   distinct(.) %>% 
#   pull(gene)
# 
# phenos = phenos %>% 
#   select(Taxa, all_of(genes_in_expression_dat))
# 
# print(paste0("There..are...", length(genes_in_expression_dat), "...Phenolic..and..TF..genes...","out..of...",nrow(my_genes),"...present..in..this..data:","...",my_files))
# print(paste0("Which..is..",round(length(genes_in_expression_dat)/nrow(my_genes),2),"%"))

print(paste0("There..are...",ncol(phenos)-1,"...genes present..in..this..data..after..filtering:","...",my_files))

## Read in snp PCs and PEER factors and bind an NA column for fitting models with 0 PCs
snp_pcs = fread(paste0(my_phenos_path,my_files,"_snp_expression_covariates.txt"), data.table = FALSE) %>% 
  cbind(NA,.)

## Read in 1000 permutation thresholds for cis associations
cis_threshold = fread(paste0(my_permutations_path,my_files,".cis_qtl.txt.gz")) %>% 
  select(phenotype_id,pval_nominal_threshold)

# ## Read in permutation thresholds for trans associations
trans_threshold = fread(paste0(my_permutations_path,my_files,"_trans_qtl.txt"), header = TRUE) %>%
  setNames(., c("perm","pval"))
trans_threshold = sort(trans_threshold$pval)[ceiling(nrow(trans_threshold)*0.05)]

## Set where the results will go
setwd(paste0(my_results_path,my_files))

print("Splitting expression data into chunks")
Taxa = phenos$Taxa
phenos = phenos %>% select(-Taxa)
my_split = 2500
nc = ncol(phenos)
phenos = lapply(split(as.list(phenos), cut(1:nc, my_split, labels = FALSE)), as.data.frame)
phenos = cbind(Taxa,phenos[[sp]])

## For testing one gene now
phenos = phenos[,c(1:2)]

## Run da loop
for (i in colnames(phenos)[2:ncol(phenos)]){
  my_results = list()
  my_filtered_results = list()
  naive_results = list()
  gc()
  
  tic()
  
  fitlist=MVP(
    phe = phenos %>% select(Taxa,all_of(i)),
    geno = myGD,
    map = map,
    priority="memory",
    method = "GLM",
    file.output = FALSE)
  
    naive_results[[i]] = bind_cols(as_tibble(fitlist$map), as_tibble(fitlist$glm.results)) %>%
    setNames(., c("SNP","CHROM","POS","REF","ALT","effect_naive_model","std_err_naive_model","naive_model_pval"))
  
  
  for(j in seq(1,7,1)){
    if (j==1) {
      snp_pcs_use = NULL
    }
    if (j==2){
      snp_pcs_use = snp_pcs[,c(2:6)]
    } 
    if (j==3){
      snp_pcs_use = snp_pcs[,c(2:11)]
    }
    if (j==4){
      snp_pcs_use = snp_pcs[,c(2:16)]
    }
    if (j==5){
      snp_pcs_use = snp_pcs[,c(2:21)]
    }
    if (j==6){
      snp_pcs_use = snp_pcs[,c(2:26)]
    }
    if (j==7){
      snp_pcs_use = snp_pcs[,c(2:31)]
    }
    
    my_threshold = cis_threshold %>% filter(phenotype_id==all_of(i))
    
    if(nrow(my_threshold)==0){
      my_threshold = max(cis_threshold$pval_nominal_threshold)
    }else{
      my_threshold = my_threshold %>% pull(pval_nominal_threshold)
    }
    
    fitlist=MVP(
      phe = phenos %>% select(Taxa,all_of(i)),
      geno = myGD,
      map = map,
      CV.MLM = snp_pcs_use,
      priority="memory",
      method = "MLM",
      K = kinship,
      vc.method = "HE",
      file.output = FALSE)
    
    if(j==1){
      my_results[[i]] = bind_cols(naive_results[[i]], as_tibble(fitlist$mlm.results[,3])) %>% 
        setNames(., c("SNP","CHROM","POS","REF","ALT","effect_naive_model",
                      "std_err_naive_model","naive_model_pval", "k_model_pval"))
    }
    if (j==2){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model1_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==3){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model2_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==4){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model3_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==5){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model4_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==6){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model5_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==7){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model6_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
  }
  
  my_filtered_results[[i]] = my_results[[i]] %>%
    select(-effect_naive_model,-std_err_naive_model) %>%
    reshape2::melt(., id.vars = c("SNP","CHROM","POS","REF", "ALT")) %>%
    group_by(variable) %>%
    filter(value <= my_threshold) %>%
    ungroup() %>%
    count(SNP) %>%
    left_join(., my_results[[i]]) %>%
    mutate(gene_id = i) %>%
    mutate(cis_permutation_threshold = my_threshold, .before = effect_naive_model) %>%
    mutate(trans_permutation_threshold = trans_threshold, .before = effect_naive_model)
  
  
  # ## For saving all results...for large snp files this doesnt make sense since most results are not significant
  # my_results[[i]] = my_results[[i]] %>%
  #   mutate(Trait = i, .before = 1) %>%
  #   mutate(Permutation_threshold = my_threshold, .before = effect_naive_model)
  
  my_filtered_results = bind_rows(my_filtered_results)
  # colnames(my_filtered_full_results)[3] = "N_models_support"
  
  fwrite(my_filtered_results, paste0(my_files,"_eQTL_results_",i,"chrom",chrom,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  toc()
  print(paste0("################################################  Done with trait:", i, "...which is...",
               grep(i, colnames(phenos)), "...out of....", length(colnames(phenos)[2:ncol(phenos)]), "   ######################################"))

}

# condor_qedit <chtc-user> RequestMemory 3000
# condor_qedit <chtc-user> RequestDisk 7500
# condor_q <chtc-user> -af RequestMemory MemoryUsage
# condor_release <chtc-user>

