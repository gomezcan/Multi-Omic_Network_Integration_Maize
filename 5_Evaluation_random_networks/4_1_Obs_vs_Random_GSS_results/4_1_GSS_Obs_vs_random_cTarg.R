suppressMessages(library(parallel))
suppressMessages(library(data.table))
suppressMessages(library(scales))
suppressMessages(library(tidyverse))
suppressMessages(library(rrvgo))
suppressMessages(library(reshape2))
suppressMessages(library(GOSemSim))
suppressMessages(library(org.Zmays.eg.db))

rm(list=ls())
######################################
####           Data input         ####
######################################
# TFs of interest
AllTFs <- unique(fread("../TFs_2_test_RNets.txt", h=F)$V1)

# random obs
DFGO <- fread("Input_Random_cTarget.txt")

# Pre-calculated semantic similarity for from org.Zmays.eg.db for BP term
Zm.GOSemSim.BP <- readRDS("../Zm.GOSemSim.BP.rds")

# Read observed GOs: True values
ObsGOs <- fread("../ReduceGOterms_All_methods.txt", h=T)


######################################

gss <- function(go1, go2) {
  val <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")
  return(val)
}


#tf = "Zm00001d005578"

Getgss_ctarg <- function(tf){
  # 
  # # GOs from random nets:  Common targets
  random_ctarg <- subset(DFGO, TF==tf) # Input_Random_cTarget.txt
  
  if(nrow(random_ctarg) == 0){
    cat(paste0(' .. Salado: ', tf, " .. \n"))
  }
  
  else{
  random_ctarg$.id <- as.character(random_ctarg$.id)
  random_ctarg <- split(random_ctarg$GO.ID, random_ctarg$.id)
  # 
  # # GOs from random nets:  Common targets
  # random_cfunct <- subset(DFGO_cfunct, TF==tf)
  # random_cfunct$.id <- as.character(random_cfunct$.id)
  # random_cfunct <- split(random_cfunct$GO.ID, random_cfunct$.id)
  
  # GOs from random nets:  Common targets
  # random_nbase <- subset(DFGO, TF==tf)
  # random_nbase$.id <- as.character(random_nbase$.id)
  # random_nbase <- split(random_nbase$GO.ID, random_nbase$.id)
  
  
  # Read observed GOs: True values
  #ObsGOs_ctarg <- subset(ObsGOs, TF==tf & Method=='Com.Target')$parent
  #ObsGOs_cfunc <- subset(ObsGOs, TF==tf & Method=='Com.Function')$parent
  ObsGOs_ctarg <- subset(ObsGOs, TF==tf & Method=='Com.Target')$parent
  
  gssvals <- lapply(random_ctarg, function(x) gss(ObsGOs_ctarg, x))
  
  gssvals <- as.data.frame(unlist(gssvals)) %>% 
    rownames_to_column() %>% 
    as.data.table()
  colnames(gssvals) <- c(".id", 'GSSr')
  gssvals[,'TF'] <- tf
  gssvals[,'Method'] <- 'Com.Target'
  
  fwrite(gssvals, paste0("GSSrDB/GSSr_cTarget_", tf,'.txt'), row.names = F, quote = F, sep = '\t')
  
  return(cat(paste0(' .. Done TF: ', tf, " ..\n")))
  } 
}

lapply(AllTFs, Getgss_ctarg)