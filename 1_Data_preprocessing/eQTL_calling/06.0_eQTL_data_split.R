## Split data to be processed by chtc faster
setwd("~/eqtl_pdi_coexpr_condor/")

k = seq(1,10,1)
sp = 40
chrom = 150

## arguments from the command line:
# k = the file name suffix - should only be 1 ... 10
# sp = the chunk being processed. This will come from the text file

## Permutations run for each dataset separately using tensorQTL
library(data.table)
library(tidyverse)
library(rMVP)

## Set paths to data in rMVP format
my_phenos_path = path.expand(paste0(getwd(),"/eqtl_pdi_coexpr_",k,"/data/07_rMVP_data/"))

## Split map files
for(i in seq_along(my_phenos_path)){
  setwd(my_phenos_path[i])
  my_files = list.files(pattern = "\\.ind$") %>% gsub(".geno.ind","",.)

  ## Read in map data and assign a sudo chromosome by splitting the map into equal parts of 
  n = chrom
  nr = length(count.fields(paste0(my_phenos_path[i],my_files,".geno.map"), sep = "\t"))
  idx = rep(1:n, each=ceiling((nr-1)/n), length.out=(nr-1))
  map = fread(paste0(my_phenos_path[i],my_files,".geno.map")) %>% 
    mutate(CHROM = gsub("chr","",CHROM)) %>% 
    mutate(idx = all_of(idx)) %>% 
    mutate(gd_idx = seq(1:(nr-1)))
  
  for(j in 1:n){
    map_tmp = map %>% filter(idx == all_of(j))
    fwrite(map_tmp, paste0(my_files,".geno.sudo.chrom",j,".map",".txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}

## Read in expression data and filter for phenolic and TF genes present
for (i in seq_along(my_phenos_path)){
  setwd(my_phenos_path[i])
  my_files = list.files(pattern = "\\.ind$") %>% gsub(".geno.ind","",.)
  phenos = fread(paste0(my_phenos_path[i],my_files,".phe"))
  Taxa = phenos$Taxa
  phenos_tmp = phenos %>% select(-Taxa)
  my_split = sp
  nc = ncol(phenos_tmp)
  phenos_tmp = lapply(split(as.list(phenos_tmp), cut(1:nc, my_split, labels = FALSE)), as.data.frame)
  
  for(j in 1:sp){
    phenos_write = cbind(Taxa,phenos_tmp[[j]])
    fwrite(phenos_write, paste0(my_phenos_path[i],my_files,".chunk",j,".phe.txt"), sep = "\t", row.names = FALSE, col.names = TRUE)
  }
}

# ## Remove files if adjusting this parameter
# for (i in seq_along(my_phenos_path)){
#   setwd(my_phenos_path[i])
#  system("rm *.phe.txt")
#  system("rm *.map.txt")
# }

# ## Remove original large files
# for (i in seq_along(my_phenos_path)){
#   setwd(my_phenos_path[i])
#  system("rm *.phe")
#  system("rm *.map")
# }

