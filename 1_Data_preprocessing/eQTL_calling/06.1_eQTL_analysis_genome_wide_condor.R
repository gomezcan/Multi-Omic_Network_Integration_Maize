# =============================================================================
# EDITORIAL NOTE (added 2026-08-18 when archiving for publication; not by the
# original author). This is the copy delivered from the original analysis
# environment. It is a SINGLE-JOB DEBUG SNAPSHOT, not the version submitted to
# the cluster:
#   * commandArgs() parsing is commented out and k/sp/chrom are hardcoded to 1;
#   * the in-script splits (n = 10 map chunks, my_split = 900 genes) do not
#     match the production job grid, which was 150 marker blocks x 40 gene
#     chunks x 8 tissues (see 06.0_eQTL_data_split.R and 06.3_identify_jobs_to_rerun.R);
#   * results are written to the pilot directory (04_eQTL_results_phenolics_and_tfs)
#     whereas the compile step reads 05_eQTL_results_genome_wide.
# The MODEL DEFINITIONS below are the production ones and are what the
# manuscript Methods describe: one naive GLM (no covariates) plus seven kinship
# MLMs whose fixed effects are none / 5 SNP PCs / 5 SNP PCs + the top 5, 10, 15,
# 20, 25 PEER factors. The HTCondor submit description and job wrapper were not
# preserved. Restore the commandArgs block and the 150 x 40 grid to re-run.
# =============================================================================
## Create a copy of the directory where all files live
# cp -r eqtl_pdi_coexpr/ eqtl_pdi_coexpr_condor

## tar files before sending to chtc
# for i in {1..10}
# do
# tar -czvf eqtl_pdi_coexpr_${i}.tar.gz eqtl_pdi_coexpr_${i}/
# done

## Copy the files to squid
# scp *.tar.gz <chtc-user>@transfer.chtc.wisc.edu:/squid/<chtc-user>/

# # ## Install R packages in the R build on condor
# install.packages("devtools")
# install.packages("tictoc")
# install.packages("reshape2")
# install.packages("R.utils")
# devtools::install_github("xiaolei-lab/rMVP")
# devtools::install_github("tidyverse/tidyverse")
# devtools::install_github("Rdatatable/data.table")


setwd("~/eqtl_pdi_coexpr_condor/")

args = commandArgs(trailingOnly = TRUE)

# k = as.numeric(args[1])
# sp = as.numeric(args[2])
# chrom = as.numeric(args[3])
k = 1
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
my_annotation_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/data/06_annotation_data/"))
my_permutations_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/results/02_eQTL_permutations/"))
my_results_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/results/04_eQTL_results_phenolics_and_tfs/"))

## Create list of file names in rMVP format
setwd(my_phenos_path)
my_files = list.files(pattern = "\\.ind$") %>% gsub(".geno.ind","",.)

## Read in snp data
myGD = attach.big.matrix(paste0(my_phenos_path,my_files,".geno.desc"))

## Read in map data and assign a sudo chromosome by splitting the map into equal parts
n = 10
map = fread(paste0(my_phenos_path,my_files,".geno.map")) %>% 
  mutate(CHROM = gsub("chr","",CHROM)) %>% 
  mutate(idx = rep(1:n, each=ceiling(nrow(.)/n), length.out=nrow(.)))
  
## Keep daters for chunk being processed
map_idx = which(map$idx == chrom)
myGD = deepcopy(myGD,rows = map_idx)
map = map %>% filter(idx == all_of(chrom))
gc()

## Read in Kinship
kinship = attach.big.matrix(paste0(my_phenos_path,my_files,".kin.desc"))

## Read in expression data and filter for phenolic and TF genes present
phenos = fread(paste0(my_phenos_path,my_files,".phe"))

print(paste0("There..are...",ncol(phenos)-1,"...genes present..in..this..data..after..filtering:","...",my_files))

## Read in snp PCs and PEER factors and bind an NA column for fitting models with 0 PCs
snp_pcs = fread(paste0(my_phenos_path,my_files,"_snp_expression_covariates.txt"), data.table = FALSE) %>% 
  cbind(NA,.)

## Read in 1000 permutation thresholds for cis associations
cis_threshold = fread(paste0(my_permutations_path,my_files,".cis_qtl.txt.gz")) %>% 
  select(phenotype_id,pval_nominal_threshold)

## Read in permutation thresholds for trans associations
trans_threshold = fread(paste0(my_permutations_path,my_files,"_trans_qtl.txt"), header = TRUE) %>%
  setNames(., c("perm","pval"))
trans_threshold = sort(trans_threshold$pval)[ceiling(nrow(trans_threshold)*0.05)]

## Read in gene location file
gene_anno = fread(paste0(my_annotation_path,"Zea_mays.B73_RefGen_v4.50.gff3.txt"))

## Set where the results will go
setwd(paste0(my_results_path,my_files))

print("Splitting expression data into chunks")
Taxa = phenos$Taxa
phenos = phenos %>% select(-Taxa)
my_split = 900
nc = ncol(phenos)
phenos = lapply(split(as.list(phenos), cut(1:nc, my_split, labels = FALSE)), as.data.frame)
phenos = cbind(Taxa,phenos[[sp]])

## For testing one gene now
# phenos = phenos[,c(1:2)]

## Run da loop
for (i in colnames(phenos)[2:ncol(phenos)]){
  my_results = list()
  naive_results = list()
  gc()
  
  tic()
  
  fitlist=MVP(
    phe = phenos %>% select(Taxa,all_of(i)),
    geno = myGD,
    map = map %>% select(-idx),
    priority="memory",
    method = "GLM",
    file.output = FALSE)
  
    naive_results[[i]] = bind_cols(as_tibble(fitlist$map), as_tibble(fitlist$glm.results)) %>%
    setNames(., c("snp","snp_chrom","snp_pos","ref_allele","alt_allele","effect_naive_model",
                  "std_err_naive_model","naive_model_pval"))
  
  
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
        setNames(., c("snp","snp_chrom","snp_pos","ref_allele","alt_allele","effect_naive_model",
                      "std_err_naive_model","naive_model_pval", "k_model_pval"))
    }
    if (j==2){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_model_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==3){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model1_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==4){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model2_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==5){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model3_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==6){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model4_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
    if (j==7){
      pvalue = as_tibble(fitlist$mlm.results[,3]) %>% setNames("k_q_p_model5_pval")
      my_results[[i]] = bind_cols(my_results[[i]],pvalue)
    }
  }
  
  my_cis_results = my_results[[i]] %>%
    mutate(gene_id = i) %>%
    left_join(., gene_anno) %>% 
    filter(snp_chrom==gene_chrom) %>% 
    filter(abs(snp_pos-gene_TSS)<=50e3) %>% 
    select(-effect_naive_model,-std_err_naive_model) %>%
    reshape2::melt(., id.vars = c("snp","snp_chrom","snp_pos","ref_allele", "alt_allele","gene_id","gene_chrom","gene_TSS","gene_TTS")) %>%
    group_by(variable) %>%
    filter(value <= my_threshold) %>%
    ungroup() %>%
    count(snp) %>%
    left_join(., my_results[[i]]) %>%
    mutate(gene_id = i) %>%
    left_join(., gene_anno) %>% 
    mutate(cis_permutation_threshold = my_threshold, .before = effect_naive_model) %>%
    mutate(trans_permutation_threshold = trans_threshold, .before = effect_naive_model) %>% 
    mutate(snp_distance_to_TSS = abs(snp_pos-gene_TSS)) %>% 
    mutate(snp_up_down_stream = ifelse(snp_pos <= gene_TSS,paste0("upstream"),
                                         ifelse(snp_pos >= gene_TSS & snp_pos <= gene_TTS,"within_gene_body",
                                                ifelse(snp_pos >= gene_TSS, paste0("downstream"), no = NA)))) %>% 
    select(snp, snp_chrom,snp_pos, ref_allele, alt_allele, n, cis_permutation_threshold, trans_permutation_threshold, effect_naive_model,
           std_err_naive_model, naive_model_pval, k_model_pval, k_q_model_pval, k_q_p_model1_pval, k_q_p_model2_pval, k_q_p_model3_pval, k_q_p_model4_pval,
           k_q_p_model5_pval, gene_id, gene_chrom, gene_TSS, gene_TTS,snp_distance_to_TSS,snp_up_down_stream)
  colnames(my_cis_results)[6] = "n_models_support"
    
  
  my_trans_results = my_results[[i]] %>%
    mutate(gene_id = i) %>%
    left_join(., gene_anno) %>% 
    # filter(snp_chrom==gene_chrom) %>% 
    # filter(abs(snp_pos-gene_TSS)<=50e3) %>% 
    select(-effect_naive_model,-std_err_naive_model) %>%
    reshape2::melt(., id.vars = c("snp","snp_chrom","snp_pos","ref_allele", "alt_allele","gene_id","gene_chrom","gene_TSS","gene_TTS")) %>%
    group_by(variable) %>%
    filter(value <= trans_threshold) %>%
    ungroup() %>%
    count(snp) %>%
    left_join(., my_results[[i]]) %>%
    mutate(gene_id = i) %>%
    left_join(., gene_anno) %>% 
    mutate(cis_permutation_threshold = my_threshold, .before = effect_naive_model) %>%
    mutate(trans_permutation_threshold = trans_threshold, .before = effect_naive_model) %>% 
    mutate(snp_distance_to_TSS = ifelse(snp_chrom == gene_chrom, abs(snp_pos-gene_TSS),no = NA)) %>%
    mutate(snp_up_down_stream = ifelse(snp_pos <= gene_TSS,paste0("upstream"),
                                       ifelse(snp_pos >= gene_TSS & snp_pos <= gene_TTS,"within_gene_body",
                                              ifelse(snp_pos >= gene_TSS, paste0("downstream"),no = NA)))) %>% 
    select(snp, snp_chrom,snp_pos, ref_allele, alt_allele, n, cis_permutation_threshold, trans_permutation_threshold, effect_naive_model,
           std_err_naive_model, k_model_pval, k_q_model_pval, k_q_p_model1_pval, k_q_p_model2_pval, k_q_p_model3_pval, k_q_p_model4_pval,
           k_q_p_model5_pval, gene_id, gene_chrom, gene_TSS, gene_TTS,snp_distance_to_TSS,snp_up_down_stream)
  colnames(my_trans_results)[6] = "n_models_support"
  
  fwrite(my_cis_results, paste0("cis_",my_files,"_eQTL_results_",i,"chrom",chrom,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  fwrite(my_trans_results, paste0("trans_",my_files,"_eQTL_results_",i,"chrom",chrom,".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  toc()
  print(paste0("################################################  Done with trait:", i, "...which is...",
               grep(i, colnames(phenos)), "...out of....", length(colnames(phenos)[2:ncol(phenos)]), "   ######################################"))

}

# condor_qedit <chtc-user> RequestMemory 3000
# condor_qedit <chtc-user> RequestDisk 7500
# condor_q <chtc-user> -af RequestMemory MemoryUsage
# condor_release <chtc-user>

