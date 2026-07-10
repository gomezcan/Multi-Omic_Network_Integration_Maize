suppressMessages(library(tidyverse))
suppressMessages(library(ggrepel))
suppressMessages(library(ggpubr))
suppressMessages(library(data.table))
suppressMessages(library(reshape2))
suppressMessages(library(circlize))
suppressMessages(library(data.table))
suppressMessages(library(scales))
suppressMessages(library(purrr))
suppressMessages(library(parallel))
#
ReplaceName <- function(ids){
  # for (i in 1:nrow(Top45)){
  #   ids <- gsub(Top45$V1[i], Top45$V2[i], ids)
  # }
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

teQTL_RandomTarg <- function(filetf){
  # Read files with random networks and select targets 
  # associated with Source genes annotated as TFs
  ## random net
  rnet <- as_tibble(fread(paste("teQTL_NetsV2/", filetf, sep = "")))
  
  # reduce to source genes annoated as TF
  rnet <- subset(rnet, Source	%in% All_TFs)
  
  return(rnet)
}


add_PCC.GAN_Random <- function(tf){
  
  ## Read names of all PCC files available: total CoExp Nets
  pcc_files <-  list.files(path = "wPCC_net_only_TFs/", pattern = paste0("*",tf,'.txt'))
  
  # PCC Nets name
  PCC_Nets <- sapply(strsplit(pcc_files, split='.', fixed=TRUE), `[`, 2) # nets
  
  # Sub-sampling random net selecting only TF's random targets 
  #net <- lapply(R_GAN_Targ, subset, Source==tf) 
  net <- lapply(R_GAN_Targ, function(x)x[x$Source==tf,]) 
  
  # rename columns in random nets 
  net <- lapply(net, setNames, c("Source","GeneID")) 
  
  # Read each Co-expression files and add wPCC info to random targets
  c=1
  print(paste0(" ... Doing TF: ", tf, " .."))
  for (f in pcc_files){
    
    PCCf <- as_tibble(fread(paste("wPCC_net_only_TFs/",f, sep = "")))[,2:3]
    colnames(PCCf)[2] <- PCC_Nets[c]
    
    # apply left_join to add pcc values by random net
    net <- lapply(net, left_join, PCCf, by="GeneID") 
    
    c=c+1
  }
  # Get average wPCC amoung random targets by random network
  
  testR_m <- lapply(net, Get_PCCmeanByNet)
  testR_m <- lapply(testR_m, as.data.frame)
  testR_m <- lapply(testR_m, function(x){x[,"Net"] <- rownames(x); x})
  
  testR_m <- as_tibble(rbindlist(testR_m, idcol = T))
  colnames(testR_m) <- c("Rsample", "wPCCm", "Net")
  testR_m$wPCCm <- round(testR_m$wPCCm, 4)
  
  write.table(testR_m, paste0("GAN_R_CoExpData/MeanVals/wPCCm_GAN_Random_", tf, ".txt"), sep = "\t", quote = F, row.names = F)
  
  # Combine all wPCCs by random net into a only DF by TF to save
  net <- as_tibble(rbindlist(net, idcol = T))
  write.table(net, paste0("GAN_R_CoExpData/RanTarg_wPCC_by_TF/Targ_wPCC_GAN_Random_", tf, ".txt"), sep = "\t", quote = F, row.names = F)
  print(paste0(" ... Done TF: ", tf, " .."))
}

Get_PCCmeanByNet <- function(df){
  # wPCC mean by network in DF with targetrs as rows and wPCC-net by column 
  apply(df[,-c(1:2)], 2, mean, na.rm=TRUE)
}

##################################################
##########        Annotations       ##############
##################################################

## syntenic genes 
Syntenic <- as_tibble(read.table("Data/Annotations/Zm.v4.synteny.genes.txt", h=T, stringsAsFactors = F, sep = '\t'))
Syntenic <- subset(Syntenic, gene_synteny=='syntenic')$gene_id
# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

# TF and CoReg
TF_CoR <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))

# CornCYC
CornCYC <- as_tibble(read.table("Data/Annotations/corn_pathways.0210325.reduced.txt", h=T, stringsAsFactors = F))
CornCYC <- CornCYC[!(CornCYC$Gene.id %in% "unknown"),]
CornCYC <- subset(CornCYC, Gene.id %in% Syntenic)
CornCYC[,"Class"] <- "Enzyme"
colnames(CornCYC)[3] <- "GeneID"
CornCYC$Pathway.id <- gsub("-", "_", CornCYC$Pathway.id)
#
CornCYC$Pathway.name <- gsub("</i>", "", gsub("<i>", "", CornCYC$Pathway.name))

CornCYC_size <- as_tibble(as.data.frame(table(unique(CornCYC[,c(1,3)])$Pathway.id), stringsAsFactors = F))
colnames(CornCYC_size) <- c("PWY", "PWYSize")
##################################################

###################################################################################################
####################################            Networks          #################################
###################################################################################################

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"
PDI <- unique(PDI$Source)

# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
#CoExp <- unique(CoExp[,2:3])
CoExp <- unique(CoExp$Source)

All_TFs <- unique(c(PDI, CoExp, TF_CoR$GeneID))

# teQTL
teQTL <- as_tibble(unique(fread("../Fig_transeQTL/teQTL_NetworkFinal.10_11_2021.txt")))
colnames(teQTL)[1] <- "Source"

teQTLtf <- subset(teQTL, Source %in% All_TFs)

###################################################################################################


###################################################################### 
########                Add PCC to GAN network                ######## 
######################################################################

TF_in_wPCCDB <- as.character(read.table("wPCC_TFs_Done.txt", h=F)$V1)
TF_in_wPCCDB <- as_tibble(as.data.frame(table(TF_in_wPCCDB), stringsAsFactors = F)) 

TF_in_wPCCDB <- TF_in_wPCCDB[order(-TF_in_wPCCDB$Freq),]



# Read random GAN nets
R_GAN_files <-  list.files(path = "teQTL_NetsV2/", pattern = '*.txt')
R_GAN_Targ <- lapply(R_GAN_files, teQTL_RandomTarg)
names(R_GAN_Targ) <- gsub(".txt", "", gsub("andom_teQTL.", "", R_GAN_files))


Lgenes <- length(TF_in_wPCCDB$TF_in_wPCCDB)
#Lgenes <- 3
w=50 # Size of range to test
print(".. Ready to start ..")

for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- TF_in_wPCCDB$TF_in_wPCCDB[Start:end]
    print(listtotest)
    #print(end)
    mclapply(listtotest, add_PCC.GAN_Random, mc.cores=w)
    
  }
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- TF_in_wPCCDB$TF_in_wPCCDB[Start:max]
    print(listtotest)
    #print(max)
    mclapply(listtotest, add_PCC.GAN_Random, mc.cores=w)
  }
}

