library(data.table)
library(tidyverse)
library(rMVP)
library(RColorBrewer)

## Read in inbred line subpop information from Mazaheri et al., 2019
my_phenos_path_tmp = path.expand("~/Stalk_anatomical_and_compostion/data/03_raw_genotypic_data/genos_in_WiDiv942_SNP_text_files.csv")
my_phenos_path = path.expand("~/Stalk_anatomical_and_compostion/data/03_raw_genotypic_data/widiv_942_genoNames.csv")
subpop_tmp = fread(my_phenos_path_tmp) %>% 
  select(genotype_hapmap, genotype_mazaheri) %>% 
  filter(complete.cases(.))
subpop = fread(my_phenos_path) %>% 
  select(Genotype_Mazehari, Subpopulation) %>% 
  left_join(.,subpop_tmp, by = c("Genotype_Mazehari"="genotype_mazaheri")) %>% 
  select(genotype_hapmap, Subpopulation) %>% 
  setNames(., c("Taxa","Subpopulation")) %>% 
  filter(complete.cases(.))

## Create covariates file to use for each eQTL analysis SNP PCs and PEER factors and make plots of these covariates for each analysis
setwd("~/eqtl_pdi_coexpr/data/07_rMVP_data/")
my_pcs = list.files(pattern = "\\.pc.desc$")
my_inds = list.files(pattern = "\\.ind$")

## run da loop 
for (i in seq_along(my_inds)){
  
## SNP PCs
setwd("~/eqtl_pdi_coexpr/data/07_rMVP_data/")
pcs = attach.big.matrix(my_pcs[i])
pcs = pcs[,]
inds = fread(my_inds[i], header = FALSE) %>% 
  left_join(., subpop, by = c("V1"="Taxa")) %>% 
  mutate(Subpopulation = as.factor(Subpopulation)) %>% 
  data.frame(.)

my_prefix = gsub(".geno.ind","", my_inds[i])

## Expression PEER factors
setwd("~/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/")
peer = fread(paste0(my_prefix, ".PEER_covariates.txt"))
peer_genos = peer[,1] %>% pull(ID)
peer = t(peer[,-1])

## Put them in the same order
print(identical(as.character(inds$V1), rownames(peer)))
peer = peer[match(inds$V1,rownames(peer)),]
print(identical(as.character(inds$V1), rownames(peer)))

## Scatterplots of top 5 PCs and PEER facotrs
lbls = paste("PC", 1:5, "\n", sep="")
lbls1 = paste("PEER", 1:5, "\n", sep="")
n = 12
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
col_vector = unlist(mapply(brewer.pal, 12, "Paired"))
set.seed(1234)
pie(rep(1,n), col=sample(col_vector, n))
my_cols=sample(col_vector, n)

my_file_name = gsub(".geno.ind","", my_inds[i])

setwd("~/eqtl_pdi_coexpr/results/03_eQTL_covariates_plots/")
png(paste0(my_file_name,"_snp_covar.png"), width = 8, height = 6,  units = "in", res = 400)
par(bg="transparent")
pairs(pcs, col=my_cols[inds$Subpopulation], labels=lbls, pch = 19,  cex = 0.6,lower.panel = NULL)
dev.off()

png(paste0(my_file_name,"_expression_peer.png"), width = 8, height = 6,  units = "in", res = 400)
par(bg="transparent")
pairs(peer[,1:5], col=my_cols[inds$Subpopulation], labels=lbls1, pch = 19,  cex = 0.2,lower.panel = NULL,cex.labels = 0.55)
dev.off()

## Combine the expression and snp covariates
snp_expr_covar = cbind(pcs,peer[,c(1:25)])

## Write the file
setwd("~/eqtl_pdi_coexpr/data/07_rMVP_data/")
fwrite(snp_expr_covar, file = paste0(my_file_name,"_snp_expression_covariates.txt"), col.names = FALSE, row.names = FALSE, sep = "\t") 
}

## Make legend for plots
setwd("~/eqtl_pdi_coexpr/results/03_eQTL_covariates_plots/")
png(paste0("00_legend.png"), width = 8, height = 6,  units = "in", res = 400)
par(bg="transparent")
plot(NULL)
legend(0, 1, as.vector(unique(subpop$Subpopulation)), fill=my_cols)
dev.off()

## Make covariates dataframes for transeQTL permutation tests to be run against
for (i in seq_along(my_inds)){
  my_prefix = gsub(".geno.ind","", my_inds[i])
  
  ## Expression PEER factors
  setwd("~/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/")
  peer = fread(paste0(my_prefix, ".PEER_covariates.txt"))
  ## Keep top 25 for each dataset and write it out
  peer_reduced= peer[1:25,]
  fwrite(peer_reduced, file = paste0(my_prefix, ".PEER_covariates_top25.txt"), col.names = TRUE, row.names = FALSE, sep = "\t") 
  print(paste0("Done with...",my_prefix))
}

for (i in seq_along(my_inds)){
  my_prefix = gsub(".geno.ind","", my_inds[i])
  
  ## Expression PEER factors
  setwd("~/eqtl_pdi_coexpr/data/03_expression_data/03_tensorqtl_data/")
  peer = fread(paste0(my_prefix, ".PEER_covariates.txt"))
  ## Keep top 25 for each dataset and write it out
  peer_reduced = peer[1,]
  peer_reduced[1,] = 0
  fwrite(peer_reduced, file = paste0(my_prefix, ".PEER_covariates_null.txt"), col.names = TRUE, row.names = FALSE, sep = "\t") 
  print(paste0("Done with...",my_prefix))
}

