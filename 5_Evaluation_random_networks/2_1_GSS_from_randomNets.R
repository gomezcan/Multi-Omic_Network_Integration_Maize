suppressMessages(library(GOSemSim))
suppressMessages(library(reshape2))
suppressMessages(library(GeneOverlap))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
#library(org.Zmays.eg.db)

args = commandArgs(trailingOnly=TRUE)

# Pre-calculate semantic similarity
#Zm.GOSemSim.BP <- GOSemSim::godata(org.Zmays.eg.db, keytype="GENENAME", ont = 'BP')
#saveRDS(Zm.GOSemSim.BP, file = "Zm.GOSemSim.BP.rds") 



#####################################################
#############      Functions            #############
#####################################################
Get_SS_GRN_CEN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GRN:CEN")
  
  c= 1
  #tf=AllTFs[2]
  
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GRN, TF == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_CEN, TF == tf)$GO.ID
    #
    if (length(go1) > 0 & length(go2) > 0){
      # Test GsS only if both nets have at least a GO term
      ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")
      GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
      
      print(' .. Step 1 ..')
      # Matrix of GOs
      M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
      M_ss <- data.frame(reshape2::melt(M_ss), stringsAsFactors = F) 
      colnames(M_ss) <- c("GO_GRN", "GO_CEN", "GSS") 
      print(' .. Step 2 ..')
      
      BP_GRN_test <- subset(GOsDB_GRN, TF==tf & GO.ID %in% M_ss$GO_GRN)
      BP_CEN_test <- subset(GOsDB_CEN, TF==tf & GO.ID %in% M_ss$GO_CEN)
      
      M_ss <- left_join(M_ss, BP_GRN_test[,c(1,2,4,9)], by=c("GO_GRN"="GO.ID"))
      M_ss <- left_join(M_ss, BP_CEN_test[,c(1,2,4,9)], by=c("GO_CEN"="GO.ID"))
      colnames(M_ss) <- c("GO.GRN", "GO.CEN", "GSS", 
                          "Term.GRN",  "Sig.Terms.GRN", "FDR.GRN",
                          "Term.CEN", "Sig.Terms.CEN", "FDR.CEN")
      print(' .. Step 3 ..')
      
      namef <- paste0("BP_GSS_CommonFunction_random/GSS.",netID,'.',tf,".GRN_CEN.txt")
      cat(namef)
      cat('\n')
      write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
      print(paste0(" .. Done ", c, " TFs .."))
      c = c+1
      
    }
    else{
      GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- 0
    }
    
  }
  return(GOSIM_DF) 
}

Get_SS_GRN_GAN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GRN:GAN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GRN, TF == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_GAN, TF == tf)$GO.ID
    #
    if (length(go1) > 0 & length(go2) > 0){
      ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")
      GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
      
      # Matrix of GOs
      M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
      M_ss <- data.frame(reshape2::melt(M_ss), stringsAsFactors = F) 
      colnames(M_ss) <- c("GO_GRN", "GO_GAN", "GSS") 
      
      BP_GRN_test <- subset(GOsDB_GRN, TF==tf & GO.ID %in% M_ss$GO_GRN)
      BP_GAN_test <- subset(GOsDB_GAN, TF==tf & GO.ID %in% M_ss$GO_GAN)
      
      M_ss <- left_join(M_ss, BP_GRN_test[,c(1,2,4,9)], by=c("GO_GRN"="GO.ID"))
      M_ss <- left_join(M_ss, BP_GAN_test[,c(1,2,4,9)], by=c("GO_GAN"="GO.ID"))
      colnames(M_ss) <- c("GO.GRN", "GO.GAN", "GSS", 
                          "Term.GRN", "Sig.Terms.GRN", "FDR.GRN",
                          "Term.GAN", "Sig.Terms.GAN", "FDR.GAN")
      
      namef <- paste0("BP_GSS_CommonFunction_random/GSS.",netID,'.',tf,".GRN_GAN.txt")
      
      cat(namef)
      cat('\n')
      write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
      print(paste0(" .. Done ", c, " TFs .."))
      c = c+1
    }
    else{ GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- 0 }
  }
  return(GOSIM_DF) 
}

Get_SS_GRN_eGRN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GRN:eGRN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GRN, TF == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_eGRN, TF == tf)$GO.ID
    
    if (length(go1) > 0 & length(go2) > 0){
      ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")    
      GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
      
      # Matrix of GOs
      M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
      M_ss <- data.frame(reshape2::melt(M_ss), stringsAsFactors = F) 
      colnames(M_ss) <- c("GO_GRN", "GO_eGRN", "GSS")
      
      BP_GRN_test <- subset(GOsDB_GRN, TF==tf & GO.ID %in% M_ss$GO_GRN)
      BP_eGRN_test <- subset(GOsDB_eGRN, TF==tf & GO.ID %in% M_ss$GO_eGRN)
      
      M_ss <- left_join(M_ss, BP_GRN_test[,c(1,2,4,9)],  by=c("GO_GRN"="GO.ID"))
      M_ss <- left_join(M_ss, BP_eGRN_test[,c(1,2,4,9)], by=c("GO_eGRN"="GO.ID"))
      colnames(M_ss) <- c("GO.GRN", "GO.eGRN", "GSS", 
                          "Term.GRN",  "Sig.Terms.GRN", "FDR.GRN",
                          "Term.eGRN", "Sig.Terms.eGRN", "FDR.eGRN")
      
      namef <- paste0("BP_GSS_CommonFunction_random/GSS.",netID,'.',tf,".GRN_eGRN.txt")
      cat(namef)
      cat('\n')
      write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
      
      print(paste0(" .. Done ", c, " TFs .."))
      c = c+1
    }
    else {GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- 0}
    
  }
  return(GOSIM_DF)
  
}

Get_SS_CEN_GAN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="CEN:GAN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_CEN, TF == tf)$GO.ID
    
    # CoExp
    go2 <- subset(GOsDB_GAN, TF == tf)$GO.ID
    
    if (length(go1) > 0 & length(go2) > 0){
      ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")
      GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
      
      # Matrix of GOs
      M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
      M_ss <- data.frame(reshape2::melt(M_ss), stringsAsFactors = F) 
      colnames(M_ss) <- c("GO_CEN", "GO_GAN", "GSS")
      
      BP_CEN_test <- subset(GOsDB_CEN, TF==tf & GO.ID %in% M_ss$GO_GRN)
      BP_GAN_test <- subset(GOsDB_GAN, TF==tf & GO.ID %in% M_ss$GO_GAN)
      
      M_ss <- left_join(M_ss, BP_CEN_test[,c(1,2,4,9)], by=c("GO_CEN"="GO.ID"))
      M_ss <- left_join(M_ss, BP_GAN_test[,c(1,2,4,9)], by=c("GO_GAN"="GO.ID"))
      colnames(M_ss) <- c("GO.CEN", "GO.GAN", "GSS", 
                          "Term.CEN", "Sig.Terms.CEN", "FDR.CEN",
                          "Term.GAN", "Sig.Terms.GAN", "FDR.GAN")
      
      namef <- paste0("BP_GSS_CommonFunction_random/GSS.",netID,'.',tf,".CEN_GAN.txt")
      cat(namef)
      cat('\n')
      write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)  
      
      print(paste0(" .. Done ", c, " TFs .."))
      c = c+1
    }
    else{GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- 0}
    
  }
  return(GOSIM_DF)
  
}

Get_SS_eGRN_CEN <- function(tfslist){
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="eGRN:CEN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_eGRN, TF == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_CEN, TF == tf)$GO.ID
    
    if (length(go1) > 0 & length(go2) > 0){
      ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")    
      GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
      
      # Matrix of GOs
      M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
      M_ss <- data.frame(reshape2::melt(M_ss), stringsAsFactors = F) 
      colnames(M_ss) <- c("GO_eGRN", "GO_CEN", "GSS")
      
      BP_eGRN_test <- subset(GOsDB_eGRN, TF==tf & GO.ID %in% M_ss$GO_eGRN)
      BP_CEN_test <-  subset(GOsDB_CEN, TF==tf & GO.ID %in% M_ss$GO_CEN)
      
      M_ss <- left_join(M_ss, BP_eGRN_test[,c(1,2,4,9)], by=c("GO_eGRN"="GO.ID"))
      M_ss <- left_join(M_ss, BP_CEN_test[,c(1,2,4,9)], by=c("GO_CEN"="GO.ID"))
      colnames(M_ss) <- c("GO.eGRN", "GO.CEN", "GSS", 
                          "Term.eGRN", "Sig.Terms.eGRN", "FDR.eGRN",
                          "Term.CEN", "Sig.Terms.CEN", "FDR.CEN")
      
      
      namef <- paste0("BP_GSS_CommonFunction_random/GSS.",netID,'.',tf,".eGRN_CEN.txt")
      cat(namef)
      cat('\n')
      write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
      print(paste0(" .. Done ", c, " TFs .."))
      c = c+1
    }
    else{GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- 0}
    
  }
  return(GOSIM_DF)
  
}

Get_SS_GAN_eGRN <- function(tfslist){
  
  # make df to save output
  GOSIM_DF <- tibble(Source=tfslist, SS=0, Net="GAN:eGRN")
  
  c= 1
  for (tf in tfslist){
    # PDI
    go1 <- subset(GOsDB_GAN, TF == tf)$GO.ID
    # CoExp
    go2 <- subset(GOsDB_eGRN, TF == tf)$GO.ID
    
    if (length(go1) > 0 & length(go2) > 0){
      ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine="BMA")    
      GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- ss
      
      # Matrix of GOs
      M_ss <- mgoSim(go1, go2, semData=Zm.GOSemSim.BP, measure="Wang", combine=NULL)
      M_ss <- data.frame(reshape2::melt(M_ss), stringsAsFactors = F) 
      colnames(M_ss) <- c("GO_GAN", "GO_eGRN", "GSS")
      
      BP_GAN_test <- subset(GOsDB_GAN, TF==tf & GO.ID %in% M_ss$GO_GAN)
      BP_eGRN_test <- subset(GOsDB_eGRN, TF==tf & GO.ID %in% M_ss$GO_eGRN)
      
      M_ss <- left_join(M_ss, BP_GAN_test[,c(1,2,4,9)], by=c("GO_GAN"="GO.ID"))
      M_ss <- left_join(M_ss, BP_eGRN_test[,c(1,2,4,9)], by=c("GO_eGRN"="GO.ID"))
      colnames(M_ss) <- c("GO.GAN", "GO.eGRN", "GSS", 
                          "Term.GAN", "Sig.Terms.GAN", "FDR.GAN",
                          "Term.eGRN", "Sig.Terms.eGRN", "FDR.eGRN")
      
      namef <- paste0("BP_GSS_CommonFunction_random/GSS.",netID,'.',tf,".GAN_eGRN.txt")
      
      cat(namef)
      cat('\n')
      write.table(M_ss, namef, row.names = F, sep = '\t', quote = F)
      print(paste0(" .. Done ", c, " TFs .."))
      c = c+1
    }
    else{GOSIM_DF$SS[GOSIM_DF$Source==tf]  <- 0}
    
  }
  return(GOSIM_DF)
}

chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}

#####################################################

### reduce files to targeted TFs
# Read list of TFs targeted in the analysis
# List of TFs to test enrichment
AllTFs <- unique(fread("TFs_2_test_RNets.txt", h=F)$V1)

# Pre-calculated semantic similarity for from org.Zmays.eg.db for BP term
Zm.GOSemSim.BP <- readRDS("Zm.GOSemSim.BP.rds")

# Random net id
netID = args[1]
#netID = 1000

#### identify files in netID to read
GOsDB_GRN <- list.files(path = 'BP_results_targets_random/',  pattern = paste0("Random.GRN.",netID,"_*"))
GOsDB_eGRN <- list.files(path = 'BP_results_targets_random/', pattern = paste0("Random.eGRN.",netID,"_*"))
GOsDB_CEN <- list.files(path = 'BP_results_targets_random/',  pattern = paste0("Random.CEN.",netID,"_*"))
GOsDB_GAN <- list.files(path = 'BP_results_targets_random/',  pattern = paste0("Random.GAN.",netID,"_*"))

GOsDB_GRN <- GOsDB_GRN[grepl(paste0("Random.GRN.",netID,"_"), GOsDB_GRN, )]
GOsDB_eGRN <- GOsDB_eGRN[grepl(paste0("Random.eGRN.",netID,"_"), GOsDB_eGRN, )]
GOsDB_CEN <- GOsDB_CEN[grepl(paste0("Random.CEN.",netID,"_"), GOsDB_CEN, )]
GOsDB_GAN <- GOsDB_GAN[grepl(paste0("Random.GAN.",netID,"_"), GOsDB_GAN, )]

print(' .. Donde reading inputs ..')

tf_inGRN <- chop(GOsDB_GRN, '[_]', 2)
tf_ineGRN <- chop(GOsDB_eGRN, '[_]', 2)
tf_inCEN <- chop(GOsDB_CEN, '[_]', 2)
tf_inGAN <- chop(GOsDB_GAN, '[_]', 2)

# keep tfs in all four nets
tested_TF <- subset(as.data.table(table(c(tf_inGRN, tf_ineGRN, tf_inGAN, tf_inCEN))), N==4)$V1

# remove netid from tf names
AllTFs <- AllTFs[AllTFs %in% tested_TF]

# Keep TFs 
GOsDB_GRN <- unlist(lapply(AllTFs, function(x) GOsDB_GRN[grepl(x, GOsDB_GRN)]))
GOsDB_eGRN <- unlist(lapply(AllTFs, function(x) GOsDB_eGRN[grepl(x, GOsDB_eGRN)]))
GOsDB_GAN <- unlist(lapply(AllTFs, function(x) GOsDB_GAN[grepl(x, GOsDB_GAN)]))
GOsDB_CEN <- unlist(lapply(AllTFs, function(x) GOsDB_CEN[grepl(x, GOsDB_CEN)]))

# keep tfs in all four nets
tested_TF <- subset(as.data.table(table(c(tf_inGRN, tf_ineGRN, tf_inGAN, tf_inCEN))), N==4)$V1

# remove netid from tf names
tested_TF <- gsub(paste0(netID, '_'), '', tested_TF)
AllTFs <- AllTFs[AllTFs %in% tested_TF]

# Keep TFs 
GOsDB_GRN <- unlist(lapply(AllTFs, function(x) GOsDB_GRN[grepl(x, GOsDB_GRN)]))
GOsDB_eGRN <- unlist(lapply(AllTFs, function(x) GOsDB_eGRN[grepl(x, GOsDB_eGRN)]))
GOsDB_GAN <- unlist(lapply(AllTFs, function(x) GOsDB_GAN[grepl(x, GOsDB_GAN)]))
GOsDB_CEN <- unlist(lapply(AllTFs, function(x) GOsDB_CEN[grepl(x, GOsDB_CEN)]))

# read files
GOsDB_GRN <- lapply(GOsDB_GRN, function(x) fread(paste0("BP_results_targets_random/",x)))
GOsDB_eGRN <- lapply(GOsDB_eGRN, function(x) fread(paste0("BP_results_targets_random/",x)))
GOsDB_CEN <- lapply(GOsDB_CEN, function(x) fread(paste0("BP_results_targets_random/",x)))
GOsDB_GAN <- lapply(GOsDB_GAN, function(x) fread(paste0("BP_results_targets_random/",x)))


## Combine DF results 
GOsDB_GRN <- rbindlist(GOsDB_GRN, idcol = F)
GOsDB_eGRN <- rbindlist(GOsDB_eGRN, idcol = F)
GOsDB_CEN <- rbindlist(GOsDB_CEN, idcol = F)
GOsDB_GAN <- rbindlist(GOsDB_GAN, idcol = F)

print(' .. Donde reading inputs 2 ..')


##  add FDR and filter
GOsDB_GRN <- GOsDB_GRN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(TF) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

GOsDB_eGRN <- GOsDB_eGRN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(TF) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

GOsDB_CEN <- GOsDB_CEN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(TF) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

GOsDB_GAN <- GOsDB_GAN %>% 
  #filter(Significant >= 1) %>% # Remove GOs reported without targets
  group_by(TF) %>%
  mutate("FDR" = p.adjust(classic, method = 'fdr')) %>%
  filter(FDR <= 0.1) 

print(' .. Donde FDR filters ..')

SS_GRN_CEN  <- Get_SS_GRN_CEN(AllTFs[1:3])
SS_GRN_GAN  <- Get_SS_GRN_GAN(AllTFs[1:3])
SS_GRN_eGRN <- Get_SS_GRN_eGRN(AllTFs[1:3])
SS_CEN_GAN  <- Get_SS_CEN_GAN(AllTFs[1:3])
SS_eGRN_CEN <- Get_SS_eGRN_CEN(AllTFs[1:3])
SS_GAN_eGRN <- Get_SS_GAN_eGRN(AllTFs[1:3])

print(' .. Donde GSS calculation ..')

## Combine all results from netID and save results
randomGSS <- rbind(SS_GRN_CEN,
                   SS_GRN_GAN,
                   SS_GRN_eGRN,
                   SS_CEN_GAN,
                   SS_eGRN_CEN,
                   SS_GAN_eGRN)

write.table(randomGSS, paste0('GSS_targets_random/Random_GSS.', netID, '.txt'), 
            row.names = F, col.names = T, sep = '\t',
            quote = F)