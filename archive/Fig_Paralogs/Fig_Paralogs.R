###################################################################################
#######                             Libraries                               #######
###################################################################################
library(patchwork)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(circlize)
library(factoextra)
library(reshape2)
library(GeneOverlap)
library(fgsea)
library(parallel)
library(hrbrthemes)
set.seed(42)

###################################################################################
#######                        Functions                                    #######
###################################################################################

ReplaceName <- function(ids){
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    #ids <- gsub("_", " ", ids)
  }
  return(ids)
}
#
list_to_DF <- function(list){
  # This function get a list of GOs and GeneIds to produce a table
  Net <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 2)))
  colnames(Net) <- c("GO.ID", "GeneID")
  
  GOs <- names(list)
  
  for (n in GOs){
    tem <- tibble(GO.ID=n, GeneID=list[[n]])
    Net <- rbind(Net, tem)
  }
  Net <- unique(Net)
  return(Net)
}

#
add.PCC.Info <- function(tf){
  
  # PDI DF
  net <- unique(subset(full.net, TF==tf)[,1:2])
  net[,1:2] <- apply(net[,1:2], 2, as.character)
  
  # define targets of TF
  list.targets <- unique(as.character(net$Target))
  
  # print(net)
  #### Cor PCCw
  # tf="Csa01g013850"
  file_all_pcc <-  paste("DB_wPCC/wPCC.All.",tf,".txt.gz", sep = "")
  file_seed_pcc <-  paste("DB_wPCC/wPCC.Seed.",tf,".txt.gz", sep = "")
  file_Notseed_pcc <-paste("DB_wPCC/wPCC.NotSeed.",tf,".txt.gz", sep = "")
  
  # Empty vector to save list of targets co-expressed 
  targ.coexp <- c()
  
  # Read co-expresion files
  
  # if (file.exists(file_all_pcc)==TRUE) {
  #   CorAll <- read.table(file_all_pcc, h=T)
  #   CorAll$GeneID <- as.character(CorAll$GeneID)
  #   colnames(CorAll)[1] <- "Target"
  #   
  #   ## get tail values from co-expression distribution
  #   q95 <- as.numeric(quantile(CorAll$wPCC.All, 0.95)) # top positive co-expressed
  #   q5  <- as.numeric(quantile(CorAll$wPCC.All, 0.05)) # top negative co-expressed
  #   #
  #   top.genes <- subset(CorAll, wPCC.All <= q5 | wPCC.All >= q95)
  #   list.targets.coexp <- as.character(subset(top.genes, Target %in% list.targets)$Target)
  #   print(length(list.targets.coexp))
  #   # Save co-expressed targets
  #   targ.coexp <- c(targ.coexp, list.targets.coexp)
  #   #
  #   net <- left_join(net, CorAll[,c(1:2)], by="Target")
  # }
  # 
  if (file.exists(file_seed_pcc)==TRUE) {
    #
    CorSeed <- read.table(file_seed_pcc, h=T)
    colnames(CorSeed)[1:2] <- c("Target", "wPCC.seed")
    CorSeed$Target <- as.character(CorSeed$Target)
    
    ## get tail values from co-expression distribution
    q95 <- as.numeric(quantile(CorSeed$wPCC.seed, 0.95)) # top positive co-expressed
    q5  <- as.numeric(quantile(CorSeed$wPCC.seed, 0.05)) # top negative co-expressed
    # Top genes in seed data
    top.genes <- subset(CorSeed, wPCC.seed <= q5 | wPCC.seed >= q95)
    list.targets.coexp <- as.character(subset(top.genes, Target %in% list.targets)$Target)
    print(length(list.targets.coexp))
    # Save co-expressed targets
    targ.coexp <- c(targ.coexp, list.targets.coexp)
    
    
    net <- left_join(net, CorSeed[,c(1:2)], by="Target")
  }
  
  if (file.exists(file_Notseed_pcc)==TRUE) {
    #
    CorNotSeed <- read.table(file_Notseed_pcc, h=T)
    colnames(CorNotSeed)[1:2] <- c("Target", "wPCC.NotSeed")
    CorNotSeed$Target <- as.character(CorNotSeed$Target)
    
    ## get tail values from co-expression distribution
    q95 <- as.numeric(quantile(CorNotSeed$wPCC.NotSeed, 0.95)) # top positive co-expressed
    q5  <- as.numeric(quantile(CorNotSeed$wPCC.NotSeed, 0.05)) # top negative co-expressed
    # Top genes in seed data
    top.genes <- subset(CorNotSeed, wPCC.NotSeed <= q5 | wPCC.NotSeed >= q95)
    list.targets.coexp <- as.character(subset(top.genes, Target %in% list.targets)$Target)
    print(length(list.targets.coexp))
    # Save co-expressed targets
    targ.coexp <- c(targ.coexp, list.targets.coexp)
    #
    net <- left_join(net, CorNotSeed[,c(1:2)], by="Target")
  }
  
  
  # subset only targets co-expressed
  net <- subset(net, Target %in% targ.coexp)
  print("... Done reading .. ")
  return(net)
}

PCC.file.vector <- function(file){
  
  PCC <- as_tibble(fread(file))[,2:3]
  colnames(PCC)[1:2] <- c("Target", "wPCC")
  
  # set Rank input
  Rank <- PCC$wPCC
  names(Rank) <- PCC$Target
  return(Rank)
}

#
PCC.fgsea <- function(tf){
  
  # PDI net
  net <- unique(subset(full.net, TF==tf)[,1:2])
  net[,1:2] <- apply(net[,1:2], 2, as.character)
  
  # combine Targets from GRN and CEN networks
  Targets <- unique(as.character(net$Target))
  
  #### GO file with total targets by GO
  GOs_GRN <-  paste("GOs_DB/BP_PDI.Genes.",tf,".txt", sep = "")
  GOs_CEN <-  paste("GOs_DB/BP_CoExp.Genes.",tf,".txt", sep = "")
  
  GOs_GRN <- unique(as_tibble(fread(GOs_GRN))[,1:2])
  GOs_CEN <- unique(as_tibble(fread(GOs_CEN))[,1:2])
  
  # GOs terms to be testted by TF 
  GO_targ <- subset(GOs_DB_DF, GeneID == tf)$GO
  
  # Reduce to targets if interest
  GOs_GRN <- subset(GOs_GRN, GO.ID %in% GO_targ)
  GOs_CEN <- subset(GOs_CEN, GO.ID %in% GO_targ)
  
  # combine targets and its annotations
  GOs_Net <- unique(rbind(GOs_GRN, GOs_CEN))
  #GOs_Net <- unique(GOs_GRN)
  #print(as.data.frame(table(GOs_Net$GO.ID)))
  
  # list of target genes by GO
  target.list <- split(GOs_Net$GeneID, GOs_Net$GO.ID)
  
  ## Cor wPCC file
  pcc_files <-  list.files(path = "wPCC_TF_Files/", pattern = paste0("*",tf,'.txt'))
  
  # PCC net names
  PCC_Nets <- sapply(strsplit(pcc_files, split='.', fixed=TRUE), `[`, 2) # nets
  
  # Ranks from co-expresion files
  InputRank <- lapply(paste0("wPCC_TF_Files/", pcc_files), PCC.file.vector)
  names(InputRank) <- PCC_Nets
  
  
  # GSEA by GO and by PCC network
  fgseaRes <- list()
  for (net in PCC_Nets){
    print(net)
    fgseaRes_n <- fgseaMultilevel(pathways = target.list, stats= InputRank[[net]], minSize  = 5, maxSize  = 500, eps=0)
    fgseaRes_n <- as_tibble(fgseaRes_n)
    fgseaRes_n['TF'] <- ReplaceName(tf)
    fgseaRes_n <- fgseaRes_n[order(fgseaRes_n$pval),]
    fgseaRes_n <- left_join(fgseaRes_n, GO_dic, by=c("pathway"="GO"))
    #fgseaRes_n <- subset(fgseaRes_n, padj<=0.05)
    fgseaRes[[net]] <- fgseaRes_n
    
  }
  
  fgseaRes <- as_tibble(rbindlist(fgseaRes, idcol = T))[,-c(5,6,9)]
  write.table(fgseaRes, paste0("GSEA_results/GSEA_GOs.", tf, ".txt"), quote = F, sep = "\t", row.names = F)
  
}
#
MakeExpDB <- function(tf){
  #
  exp_files <- list.files(path = "ExpDB/", pattern = "*.tsv")
  #
  out <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 2)))
  colnames(out) <- c("Net", "CPM")
  
  for (f in exp_files){
    #
    expfile = gsub(".tsv", "", f)
    # read ex file
    Exp <- as_tibble(fread(paste0("ExpDB/", f)))
    infile <- (tf %in% Exp$gid)
    
    if( infile ==TRUE){
      #
      vals <- as.numeric(Exp[Exp$gid==tf, c(2:ncol(Exp))])
      tem <- tibble(Net=expfile, CPM=vals)
      out <- rbind(out, tem)
    }
    
  }
  return(out)
}

MakeTableRankingGOs <- function(tf){
  df <- fgsea_GRN_CEN[[tf]]
  df <- subset(df, padj<=0.1)
  colnames(df)[1] <- "Network"
  sub
  df <- df %>%
    group_by(Network) %>%
    mutate(RankGO=rank(-NES)) %>%
    arrange(RankGO)
  
  GOs_DB_tem <- unique(subset(GOs_DB, .id == tf)[,c(2,3,6,8)])
  
  df <- left_join(df, unique(GOs_DB_tem[,c(1,3)]), by=c("pathway"="GO.GRN"))  
  df <- left_join(df, unique(GOs_DB_tem[,c(2,4)]), by=c("pathway"="GO.CEN"))
  
  df[,3:4] <- round(df[,3:4], 4)
  df[,5] <- round(df[,5], 2)
  
  return(df)
} 

GetZ_tissues <- function(tf) {
  
  Exp.DB <- TFs_Exp[[tf]]
  
  listTissues <- unique(Exp.DB$Net)
  
  # Define DF to save results
  dfout <- as.data.frame(matrix(0, nrow = length(listTissues), ncol = 2))
  colnames(dfout) <- c("Net", "Z")
  #
  
  IQR_all <-  IQR(Exp.DB$CPM)     # IQR
  Median_all <- mean(Exp.DB$CPM)  #
  
  # print(head(dfout))
  c=1
  for (t in listTissues){
    #
    Valmedian <- IQR(subset(Exp.DB, Net==t)$CPM) # get median along rows (genes) in Sample
    #
    dfout[c,2] <- (Valmedian - Median_all)/IQR_all
    dfout[c,1] <- t
    
    print(paste(" .. Done tissue:", t, " .. ", sep = ""))
    c=c+1
  }
  
  dfout
  
  return(dfout)
}

PCC.fgsea.plot <- function(tf, mode){
  
  # PDI net
  net <- unique(subset(full.net, TF==tf)[,1:2])
  net[,1:2] <- apply(net[,1:2], 2, as.character)
  
  # full list of genes by LRG-process
  Full.list <- split(Lipids$GeneID, Lipids$Pathway)
  
  # target list of genes by LRG-process
  Lipids_targ <- subset(Lipids, GeneID %in% unique(as.character(net$Target)))
  print(dim(Lipids_targ))
  target.list <- split(Lipids_targ$GeneID, Lipids_targ$Pathway)
  
  #### Cor PCCw file
  # tf="Csa01g013850"
  file_all_pcc <-  paste("DB_wPCC/wPCC.All.",tf,".txt.gz", sep = "")
  #file_seed_pcc <-  paste("DB_wPCC/wPCC.Seed.",tf,".txt.gz", sep = "")
  #file_Notseed_pcc <-paste("DB_wPCC/wPCC.NotSeed.",tf,".txt.gz", sep = "")
  
  # Ranks file
  InputRank <- vector()
  if (nrow(Lipids_targ)>25) {
    # Read co-expresion files
    if (file.exists(file_all_pcc)==TRUE) {
      #
      Corfile <- read.table(file_all_pcc, h=T, stringsAsFactors = F)
      colnames(Corfile)[1:2] <- c("Target", "wPCC")
      # set Rank input
      InputRank <- Corfile$wPCC
      names(InputRank) <- Corfile$Target
    }
    
    fgseaRes_full <- fgsea(pathways = Full.list, stats= InputRank, nperm = 1000, minSize  = 5, maxSize  = 900)
    fgseaRes_targ <- fgsea(pathways = target.list, stats= InputRank, nperm = 1000, minSize  = 5, maxSize  = 900)
    #
    fgseaRes <- fgseaRes_full
    fgseaRes <- fgseaRes
    
    if (mode=='targets'){
      fgseaRes <- fgseaRes_targ
      fgseaRes <- fgseaRes[order(pval), 1:7]
      
    }
    fgseaRes <- as_tibble(fgseaRes)
    
    return(list(targetlist=target.list, InputRank=InputRank, gseq=fgseaRes))
  }
  return(list(targetlist=Full.list, InputRank=InputRank, gseq=fgseaRes))
}


# Get flatten Matrix
flattenMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    Pval  =(cormat)[ut],
    Common = pmat[ut]
  )
}


##############################################################################
##################         Read data input           #########################
##############################################################################
# TFdb
TFdb <- as_tibble(read.table("Data/Annotations/TF_CoR_Mazie.txt", h=T, stringsAsFactors = F))

# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))


# Read list of TDG
TDG <- read.table("Data/Annotations/Tandem_duplicates.csv", h=F, sep = '\t') 

TDG_Table <- as_tibble(as.data.frame(matrix(0, nrow = 0, ncol = 2)))
colnames(TDG_Table) <- c("TDG_Set", "GeneID")

for (n in seq(1, nrow(TDG))) {
  g= strsplit(TDG[n,], split = ",")

  TDG_Table <- rbind(TDG_Table, 
                     tibble(TDG_Set=paste0("TDG_", n), GeneID=g[[1]])
                     )
}

TDG_Table

################################################################
########             Count TFs with a TDG               ######## 
################################################################

# unique TFS
GRN_TFs <- tibble(GeneID= unique(PDI$Source))
# add TDG infor
GRN_TFs <- left_join(GRN_TFs, TDG_Table, by=c('GeneID')) 
# selected TFs with at least 2 duplicated in dataset
GRN_TFs <- subset(GRN_TFs, TDG_Set %in% subset(as.data.frame(table(GRN_TFs$TDG_Set), stringsAsFactors = F), Freq>1)$Var1)

# unique TFS
CEN_TFs <- tibble(GeneID= unique(CoExp$Source))
# add TDG infor
CEN_TFs <- left_join(CEN_TFs, TDG_Table, by=c('GeneID')) 
# selected TFs with at least 2 duplicated in dataset
CEN_TFs <- subset(CEN_TFs, TDG_Set %in% subset(as.data.frame(table(CEN_TFs$TDG_Set), stringsAsFactors = F), Freq>1)$Var1)

# unique TFS
GAN_TFs <- tibble(GeneID= unique(teQTL$Source))
# add TDG infor
GAN_TFs <- left_join(GAN_TFs, TDG_Table, by=c('GeneID')) 
# selected TFs with at least 2 duplicated in dataset
GAN_TFs <- subset(GAN_TFs, TDG_Set %in% subset(as.data.frame(table(GAN_TFs$TDG_Set), stringsAsFactors = F), Freq>1)$Var1)

# unique TFS
PDI_TFs <- tibble(GeneID= unique(PDI$Source))
# add TDG infor
PDI_TFs <- left_join(PDI_TFs, TDG_Table, by=c('GeneID')) 
# selected TFs with at least 2 duplicated in dataset
PDI_TFs <- subset(PDI_TFs, TDG_Set %in% subset(as.data.frame(table(PDI_TFs$TDG_Set), stringsAsFactors = F), Freq>1)$Var1)

# Count common target between TFs
list_PDI_TFs <- split(PDI$Target, PDI$Source)

# Count common target between TFs
CoExp <- unique(CoExp[,2:3])
list_CEN_TFs <- split(CoExp$Target, CoExp$Source)

# Count common target between TFs
list_GAN_TFs <- split(teQTL$Target, teQTL$Source)

CountCommonTargets <- function(list.tf){
  # Count common target and test if significant
  
  go.obj <- newGOM(list.tf, list.tf, genome.size=45546) # all vs all 
  #
  Pval <- getMatrix(go.obj, name="pval")
  #J <- getMatrix(go.obj, name="Jaccard")
  common <- getMatrix(go.obj, name="intersection")
  
  out <- flattenMatrix(Pval, common)
  colnames(out)[1:2] <- c("TF1", "TF2")

  return(as_tibble(out))
}

## Get common targets and P value
CommonDF_GRN <- CountCommonTargets(list_PDI_TFs)
CommonDF_CEN <- CountCommonTargets(list_CEN_TFs)
CommonDF_GAN <- CountCommonTargets(list_GAN_TFs)

CommonDF_CEN
# Count total targets and add to DF
GRN_Total <- as_tibble(as.data.frame(table(PDI$Source)))
colnames(GRN_Total) <- c("GeneID", "Total.targ")

CEN_Total <- as_tibble(as.data.frame(table(CoExp$Source), stringsAsFactors = F))
colnames(CEN_Total) <- c("GeneID", "Total.targ")


### Add total targets ###
CommonDF_GRN <- left_join(CommonDF_GRN, GRN_Total, by= c("TF1"="GeneID"))
CommonDF_GRN <- left_join(CommonDF_GRN, GRN_Total, by= c("TF2"="GeneID"))

CommonDF_CEN <- left_join(CommonDF_CEN, CEN_Total, by= c("TF1"="GeneID"))
CommonDF_CEN <- left_join(CommonDF_CEN, CEN_Total, by= c("TF2"="GeneID"))

### Add TDG label  ###
CommonDF_GRN <- left_join(CommonDF_GRN, TDG_Table, by=c("TF1"="GeneID"))
CommonDF_GRN <- left_join(CommonDF_GRN, TDG_Table, by=c("TF2"="GeneID"))

CommonDF_CEN <- left_join(CommonDF_CEN, TDG_Table, by=c("TF1"="GeneID"))
CommonDF_CEN <- left_join(CommonDF_CEN, TDG_Table, by=c("TF2"="GeneID"))


CommonDF_GRN_TD <- subset(CommonDF_GRN, TDG_Set.x==TDG_Set.y)
CommonDF_CEN_TD <- subset(CommonDF_CEN, TDG_Set.x==TDG_Set.y)

subset(CommonDF_CEN_TD, Pval <= 0.05)
