library(hrbrthemes)
library(scales)
library(tidyverse)
library(data.table)
library(ggVennDiagram)
library(GeneOverlap)
library(topGO)
library(purrr)
library(gplots)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(viridis)
library(patchwork)
library(reshape2)
library(rrvgo)
library(org.Zmays.eg.db)
library(UpSetR)
library(fgsea)
library(GOSemSim)


################################################
##        Read enrichment results: GOs        ##
################################################
# TF names
TFdic <- as_tibble(read.table("../../Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))

GO_CommTarg <- fread("../CommonTarg_GO_enrichment.txt") %>%
  group_by(TF) %>%
  mutate(FDR=p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1)

# Common function based on GOs
GO_CommFunt <- fread("../CommonFunction_GO_enrichment.txt") %>%
  filter(FDR.1 <= 0.1 & FDR.2 <=0.1 & GSS >= 0.6)  

# Extract each pair of GOs that pass GSS filter
## Part 1
GO_CommFunt %>%
  dplyr::select(GO1, FDR.1, TF) -> GO_CommFunt1
colnames(GO_CommFunt1) <- c('GO.ID',"FDR", "TF")

## Part 2
GO_CommFunt %>%
  dplyr::select(GO2, FDR.2, TF) -> GO_CommFunt2 
colnames(GO_CommFunt2) <- c('GO.ID',"FDR", "TF")

# Combined GOs 
GO_CommFunt <- unique(rbind(GO_CommFunt1, GO_CommFunt2))

#
GO_Network <- fread("../NetworkBased_GO_Clusters_enrichment.txt") %>%
  group_by(TF) %>%
  mutate(FDR=p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1)
GO_Network

## GO enriched on DEGs
DEGs_1_GOs_TFs <- unique(DEGs_1_GOs$TF)

DEGs_1_GOs <- read.table('../DEGs_GO.db.txt', header = T)
DEGs_1_GOs_TFs <- unique(DEGs_1_GOs$TF)


# Pre-calculate semantic similarity
Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')

# Get 
set.seed(123)

chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$V2[i], TFdic$V1[i], ids)
  }
  return(ids)
}

GetGSS_DEGs_Random_CTarg <- function(tf){
  # make df to save output
  # DEGs
  go_obs <- subset(DEGs_1_GOs, TF == tf)$GO.ID
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # ComTarget
  n=length(subset(GO_CommTarg, TF == tfid)$GO.ID)
  
  print(paste0(ReplaceName(tfid), ' ', n))
  
  times=10
  
  if( n > 1 ){
    # get real GOs
    ReadGOs <- subset(GO_CommTarg, TF == tfid)$GO.ID
    
    # get random GOs
    dbCTarg <- lapply(seq(1:times), 
                      function(x) sample(GO_CommTarg$GO.ID, size = n, replace = T))
    names(dbCTarg) <- paste0('R', seq(1:times))
    
    # GSS with random GOs
    RandomGSS <- lapply(dbCTarg, function(x)
      mgoSim(go_obs, x, semData=Zm.GOSemSim.BP, 
             measure="Wang", combine='BMA'))
    
    # GSS with reads set of GOs
    RandomGSS <- as_tibble(unlist(RandomGSS))
    colnames(RandomGSS) <- c("rGSS")
    
    # real obs value
    GSS <- mgoSim(go_obs, ReadGOs, semData=Zm.GOSemSim.BP, measure="Wang", combine='BMA')
    
    # add "realGSS"
    RandomGSS[,'RealGSS'] <- GSS
    
    RandomGSS[,'Class'] <- 'Comm.Target'
    RandomGSS[,'TF'] <- tf
    return(RandomGSS) 
  }
  
  
}

GetGSS_DEGs_Random_CFunct <- function(tf){
  # make df to save output
  # DEGs
  go_obs <- subset(DEGs_1_GOs, TF == tf)$GO.ID
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # ComTarget
  n=length(subset(GO_CommFunt, TF == tfid)$GO.ID)
  
  print(paste0(ReplaceName(tfid), ' ', n))
  
  times=10
  if( n > 1 ){
    # get real GOs
    ReadGOs <- subset(GO_CommFunt, TF == tfid)$GO.ID
    
    # get random GOs
    dbCFunct <- lapply(seq(1:times), 
                       function(x) sample(GO_CommFunt$GO.ID, size = n, replace = T))
    names(dbCFunct) <- paste0('R', seq(1:times))
    
    # GSS with random GOs
    RandomGSS <- lapply(dbCFunct, function(x)
      mgoSim(go_obs, x, semData=Zm.GOSemSim.BP, 
             measure="Wang", combine='BMA'))
    
    # GSS with reads set of GOs
    RandomGSS <- as_tibble(unlist(RandomGSS))
    colnames(RandomGSS) <- c("rGSS")
    
    # real obs value
    GSS <- mgoSim(go_obs, ReadGOs, semData=Zm.GOSemSim.BP, measure="Wang", combine='BMA')
    
    # add "realGSS"
    RandomGSS[,'RealGSS'] <- GSS
    
    RandomGSS[,'Class'] <- 'Comm.Funct'
    RandomGSS[,'TF'] <- tf
    return(RandomGSS) 
  }
  
  
}

GetGSS_DEGs_Random_nbase <- function(tf){
  # make df to save output
  # DEGs
  go_obs <- subset(DEGs_1_GOs, TF == tf)$GO.ID
  
  # define tf id
  tfid <- chop(tf, '[:]', 2)
  
  # ComTarget
  n=length(subset(GO_Network, TF == tfid)$GO.ID)
  
  print(paste0(ReplaceName(tfid), ' ', n))
  
  
  if( n > 1 ){
    # get real GOs
    ReadGOs <- subset(GO_Network, TF == tfid)$GO.ID
    
    times=10
    # get random GOs
    dbnbase <- lapply(seq(1:times), 
                      function(x) sample(GO_Network$GO.ID, size = n, replace = T))
    names(dbnbase) <- paste0('R', seq(1:times))
    
    # GSS with random GOs
    RandomGSS <- lapply(dbnbase, function(x)
      mgoSim(go_obs, x, semData=Zm.GOSemSim.BP, 
             measure="Wang", combine='BMA'))
    
    # GSS with reads set of GOs
    RandomGSS <- as_tibble(unlist(RandomGSS))
    colnames(RandomGSS) <- c("rGSS")
    
    # real obs value
    GSS <- mgoSim(go_obs, ReadGOs, semData=Zm.GOSemSim.BP, measure="Wang", combine='BMA')
    
    # add "realGSS"
    RandomGSS[,'RealGSS'] <- GSS
    
    RandomGSS[,'Class'] <- 'Comm.nBase'
    RandomGSS[,'TF'] <- tf
    return(RandomGSS) 
  }
  
  
}


## Common targtes
RandomGSS_ctarg <- lapply(DEGs_1_GOs_TFs, function(x) GetGSS_DEGs_Random_CTarg(x))
mask <- unlist(lapply(RandomGSS_ctarg, function(x) is.data.frame(x)))
RandomGSS_ctarg <- rbindlist(RandomGSS_ctarg[mask])

## Common functions
RandomGSS_cFunct <- lapply(DEGs_1_GOs_TFs, function(x) GetGSS_DEGs_Random_CFunct(x))
mask <- unlist(lapply(RandomGSS_cFunct, function(x) is.data.frame(x)))
RandomGSS_cFunct <- rbindlist(RandomGSS_cFunct[mask])

## Network-based
RandomGSS_nbased <- lapply(DEGs_1_GOs_TFs, function(x) GetGSS_DEGs_Random_nbase(x))
mask <- unlist(lapply(RandomGSS_nbased, function(x) is.data.frame(x)))
RandomGSS_nbased <- rbindlist(RandomGSS_nbased[mask])

# Combined all resutls




