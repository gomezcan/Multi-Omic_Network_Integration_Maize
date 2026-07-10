###################################################################################
#######                             Libraries                               #######
###################################################################################
library(patchwork)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(scales)
library(ComplexHeatmap)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(circlize)
library(factoextra)
library(reshape2)
library(fgsea)
library(parallel)
library(hrbrthemes)
set.seed(42)

###################################################################################
#######                        Functions                                    #######
###################################################################################

ReplaceName <- function(ids){
  #   
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$TF.v4[i], TFdic$TF.Name[i], ids)
  }
  return(ids)
}

ReplaceNamePWY <- function(ids){
  
  for (i in 1:nrow(CornCYC)){
    w <- paste0('\\<', CornCYC$Pathway.id[i], '\\>')
    ids <- gsub(w, CornCYC$Pathway.name[i], ids)
    # ids <- gsub("_", " ", ids)
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

##############################################################################
##################         Read data input           #########################
##############################################################################
# TF names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=T, stringsAsFactors = F))

# PDI
PDI <- as_tibble(unique(fread("../Fig_PDI/PDI_NetworkFinal.10_11_2021.txt")[,c(2,3)]))
colnames(PDI)[1] <- "Source"
# CoExp
CoExp <- as_tibble(unique(fread("../Fig_Coexpression/CoExp_NetworkFinal.10_11_2021.txt")))
colnames(CoExp)[2] <- "Source"
# Combine networks
full.net <- unique(rbind(PDI, CoExp[,2:3]))
colnames(full.net)[1] <- "TF"
#rm(CoExp)
#rm(PDI)

# Common GOs from GSS analysis
GOs_DB <- as_tibble(fread("Sig_GOs_CEN_GRN.v2.txt"))[,-c(7,10)]

# GOs by TF
GOs_DB_DF <- unique(gather(GOs_DB[,1:3], value, key, -.id)[c(1,3)])
colnames(GOs_DB_DF) <- c("GeneID", "GO")
# GO term dic
GO_dic <- unique(tibble(GO=c(GOs_DB$GO.GRN, GOs_DB$GO.CEN), Term=c(GOs_DB$Term.GRN, GOs_DB$Term.CEN)))

# TFs to test
TF2test <- fread("TF_list.txt")$GeneID

################################################################
########         Count GOs enriched using GSEA          ######## 
################################################################

ReplaceName(TF2test[c(1,38,43)])
TF2test[c(1,38,43)]
mclapply(TF2test[c(1,38,43)], PCC.fgsea, mc.cores=3)

mclapply(TF2test[1:20], PCC.fgsea, mc.cores=20)
mclapply(TF2test[21:40], PCC.fgsea, mc.cores=20)
mclapply(TF2test[41:48], PCC.fgsea, mc.cores=8)

fgsea_GRN_CEN <- lapply(paste0("GSEA_results/GSEA_GOs.", TF2test, ".txt"), fread)
fgsea_GRN_CEN <- lapply(fgsea_GRN_CEN, as_tibble)
names(fgsea_GRN_CEN) <- TF2test

################################################################

################################################################
########         Count GOs enriched using GSEA          ######## 
################################################################

# Get expression values by CoExp network
TFs_Exp <- lapply(TF2test, MakeExpDB)
names(TFs_Exp) <- TF2test

# save tabe with GDEA result inlcuding target by network as well as 
# its corresponding 

for (tf in TF2test){
  write.table(MakeTableRankingGOs(tf), 
              paste0("GSEA_results_Ranks/GSEA_",ReplaceName(tf),"_",tf,"",".txt"), 
              row.names = F, 
              quote = F, sep = "\t")
}


MakeGSEA_Heatmap <- function(tf){
  
  df <- fgsea_GRN_CEN[[tf]]
  df <- subset(df, padj<=0.1)
  colnames(df)[1] <- "Network"
  
  # do not print GOs which does not have target in both networks
  GOs_DB_tem <- unique(subset(GOs_DB, .id == tf)[,c(2,3,6,8)])
  
  df <- left_join(df, unique(GOs_DB_tem[,c(1,3)]), by=c("pathway"="GO.GRN"))  
  df <- left_join(df, unique(GOs_DB_tem[,c(2,4)]), by=c("pathway"="GO.CEN"))
  
  df$Sig.Terms.GRN[is.na(df$Sig.Terms.GRN)] <- 0
  df$Sig.Terms.CEN[is.na(df$Sig.Terms.CEN)] <- 0
  
  
  df <- subset(df, Sig.Terms.GRN > 0 & Sig.Terms.CEN > 0)
  
  Dic <- unique(df[, c("pathway","Term")])
  
  
  
  print(paste0("... ", tf, " ..."))
  
  if( length(unique(df$Network)) >1 ){
    
    df <- df %>%
      select(Network, NES, pathway) 
    # dcast to calculate order
    # M_DF <- data.table::dcast(df[,1:3], Term ~ Network, value.var = "NES")
    M_DF <- as_tibble(as.data.frame(matrix(0, nrow = length(unique(df$pathway)), ncol = 1)))
    M_DF[,1] <- unique(df$pathway)
    colnames(M_DF)[1] <- "pathway"
    c=2
    nets <- unique(df$Network)
    print(length(nets))
    for (n in nets){
      #
      tem <- subset(df, Network==n)[,-c(1)]
      M_DF <- left_join(M_DF, tem, by="pathway")
      colnames(M_DF)[c] <- n
      print(paste0("... n: ", c-1, " ..."))
      print(dim(M_DF))
      c=c+1
    }
    
    M_DF[is.na(M_DF)] <- 0
    print(dim(M_DF))
    df <- gather(M_DF, value, key, -pathway)
    colnames(df) <- c("pathway", "Network", "NES")
    
    print(dim(df))
    df <- left_join(df, Dic, by="pathway")
    print(dim(df))
    
    
  }
  
  else{
    df <- df %>%
      select(Network, NES, Term) 
  }
  
  print(paste0(" . Ranking ", tf))
  df <- df %>%
    group_by(Network) %>%
    mutate(RankGO=rank(-NES)) %>%
    arrange(RankGO)  
  
  Rank_df <- df %>% 
    select(Term, RankGO) %>%
    group_by(Term) %>%
    summarise(RankT=median(RankGO)) %>%
    arrange(RankT)
  print(Rank_df)
  
  Top <- Rank_df$Term#[1:50]
  #Tail <- Rank_df$Term[(nrow(Rank_df)-20):nrow(Rank_df)]
  #print(Top),  #print(Tail)
  
  DFhot <- subset(df, Term %in% c(Top))
  #DFhot$Term <- factor(DFhot$Term, levels = c(Top, Tail))
  
  net_with_gsea <- unique(DFhot$Network)
  Exp <- TFs_Exp[[tf]]
  Exp <- subset(Exp, Net %in% net_with_gsea) # Print only expression from data expression
  
  Order_Exp <- Exp %>%
    group_by(Net) %>%
    mutate(Mean=mean(CPM)) %>%
    select(Net, Mean) %>%
    arrange(desc(Mean))  
  
  Order_Exp <- unique(Order_Exp)
  Order_Exp <- subset(Order_Exp, Mean >=0.1)$Net
  
  Exp$Net <- factor(Exp$Net, levels = Order_Exp)
  
  Plot1<- ggplot(Exp, aes(x=Net, y=CPM))+
    geom_jitter(size=0.2, alpha=0.3, width = 0.01)+
    theme_pubclean() +
    #scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), position = "right") +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1, size = 8),
          plot.subtitle = element_text(size=8)) +
    xlab("Network") +
    labs(subtitle = paste0(tf, " (", ReplaceName(tf), ")"))
  
  
  DFhot$NES[DFhot$NES > 3] <- 3
  DFhot$NES[DFhot$NES < -3] <- -3
  
  #########
  # M_DFhot <- reshape2::dcast(df[,1:3], Network ~ Term , value.var = "NES")
  # row.names(M_DFhot) <- M_DFhot$Network
  # M_DFhot <- t(M_DFhot[,-c(1)])
  # M_DFhot[is.na(M_DFhot)] <- 0
  # d <- dist(M_DFhot, method = "euclidean") # distance matrix
  # fit <- hclust(d, method="complete")
  # Label_order <- fit$labels[fit$order]
  # DFhot$Term <- factor(DFhot$Term, levels = c(rev(Label_order)))
  #########
  
  DFhot$Term <- factor(DFhot$Term, levels = c(rev(Rank_df$Term)))
  
  DFhot <- subset(DFhot, Network %in% Order_Exp) # heatmap if pass exp filter 
  DFhot$Network <- factor(DFhot$Network, levels = Order_Exp)
  
  Plot2 <- ggplot(DFhot, aes(x=Network, y=Term, fill=NES)) + 
    geom_tile() +
    geom_text(aes(label=RankGO), size=1) +
    scale_fill_viridis_c(direction = 1, option = "A", limits=c(-3,3)) +
    theme_bw()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 5)) 
  
  PlotFinal <- Plot1+ Plot2 +  plot_layout(heights = c(0.5, 5))
  #PlotFinal <- ggarrange(Plot1, Plot2, 
  #                      nrow = 2, align = 'v', heights = c(1,2))
  
  return(PlotFinal)
  
}

for (tf in TF2test){
  MakeGSEA_Heatmap(tf)
  ggsave(paste0("Figures_/GSEA_",ReplaceName(tf),"_",tf,"",".pdf"), width=8, height=18, units = c("in"))
}


######
# Size: 20x15
MakeGSEA_Heatmap(TF2test[6])
MakeGSEA_Heatmap(TF2test[38])
MakeGSEA_Heatmap(TF2test[43])

################################################################
########    Top candidates: GO and GSEA analyzes        ######## 
################################################################


MakeGSEA_HeatmapSubset <- function(tf, TargetTerms){
  
  df <- fgsea_GRN_CEN[[tf]]
  df <- subset(df, padj<=0.1)
  colnames(df)[1] <- "Network"
  
  GOs_DB_tem <- unique(subset(GOs_DB, .id == tf)[,c(2,3,6,8)])
  df <- left_join(df, unique(GOs_DB_tem[,c(1,3)]), by=c("pathway"="GO.GRN"))  
  df <- left_join(df, unique(GOs_DB_tem[,c(2,4)]), by=c("pathway"="GO.CEN"))
  
  df$Sig.Terms.GRN[is.na(df$Sig.Terms.GRN)] <- 0
  df$Sig.Terms.CEN[is.na(df$Sig.Terms.CEN)] <- 0
  
  
  df <- subset(df, Sig.Terms.GRN > 0 & Sig.Terms.CEN > 0)
  
  Dic <- unique(df[, c("pathway","Term")])
  
  
  
  print(paste0("... ", tf, " ..."))
  
  if( length(unique(df$Network)) >1 ){
    
    df <- df %>%
      select(Network, NES, pathway) 
    # dcast to calculate order
    # M_DF <- data.table::dcast(df[,1:3], Term ~ Network, value.var = "NES")
    M_DF <- as_tibble(as.data.frame(matrix(0, nrow = length(unique(df$pathway)), ncol = 1)))
    M_DF[,1] <- unique(df$pathway)
    colnames(M_DF)[1] <- "pathway"
    c=2
    nets <- unique(df$Network)
    print(length(nets))
    for (n in nets){
      #
      tem <- subset(df, Network==n)[,-c(1)]
      M_DF <- left_join(M_DF, tem, by="pathway")
      colnames(M_DF)[c] <- n
      print(paste0("... n: ", c-1, " ..."))
      print(dim(M_DF))
      c=c+1
    }
    
    M_DF[is.na(M_DF)] <- 0
    print(dim(M_DF))
    df <- gather(M_DF, value, key, -pathway)
    colnames(df) <- c("pathway", "Network", "NES")
    
    print(dim(df))
    df <- left_join(df, Dic, by="pathway")
    print(dim(df))
    
    
  }
  
  else{
    df <- df %>%
      select(Network, NES, Term) 
  }
  
  print(paste0(" . Ranking ", tf))
  df <- df %>%
    group_by(Network) %>%
    mutate(RankGO=rank(-NES)) %>%
    arrange(RankGO)  
  
  Rank_df <- df %>% 
    select(Term, RankGO) %>%
    group_by(Term) %>%
    summarise(RankT=median(RankGO)) %>%
    arrange(RankT)
  
  Top <- Rank_df$Term#[1:50]
  #Tail <- Rank_df$Term[(nrow(Rank_df)-20):nrow(Rank_df)]
  #print(Top),  #print(Tail)
  
  DFhot <- subset(df, Term %in% c(Top))
  #DFhot$Term <- factor(DFhot$Term, levels = c(Top, Tail))
  
  net_with_gsea <- unique(DFhot$Network)
  Exp <- TFs_Exp[[tf]]
  Exp <- subset(Exp, Net %in% net_with_gsea) # Print only expression from data expression
  
  Order_Exp <- Exp %>%
    group_by(Net) %>%
    mutate(Mean=mean(CPM)) %>%
    select(Net, Mean) %>%
    arrange(desc(Mean))  
  
  Order_Exp <- unique(Order_Exp)
  Order_Exp <- subset(Order_Exp, Mean >=0.1)$Net
  
  Exp$Net <- factor(Exp$Net, levels = Order_Exp)
  
  Plot1<- ggplot(Exp, aes(x=Net, y=CPM))+
    geom_jitter(size=0.2, alpha=0.3, width = 0.01)+
    theme_pubclean() +
    #scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), position = "right") +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1, size = 8),
          plot.subtitle = element_text(size=8)) +
    xlab("Network") +
    labs(subtitle = paste0(tf, " (", ReplaceName(tf), ")"))
  
  
  DFhot$NES[DFhot$NES > 3] <- 3
  DFhot$NES[DFhot$NES < -3] <- -3
  
  DFhot$Term <- factor(DFhot$Term, levels = c(rev(Rank_df$Term)))
  DFhot <- subset(DFhot, Network %in% Order_Exp) # heatmap if pass exp filter 
  DFhot$Network <- factor(DFhot$Network, levels = Order_Exp)
  print(DFhot)
  
  DFhot <- filter(DFhot, grepl(paste(TargetTerms, collapse="|"), Term)) 
  
  Plot2 <- ggplot(DFhot, aes(x=Network, y=Term, fill=NES)) + 
    geom_tile() +
    geom_text(aes(label=RankGO), size=3) +
    scale_fill_viridis_c(direction = 1, option = "A", limits=c(-3,3)) +
    theme_bw()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 7)) 
  
  PlotFinal <- Plot1+ Plot2 +  plot_layout(heights = c(1, 5))
  #PlotFinal <- ggarrange(Plot1, Plot2, 
  #                      nrow = 2, align = 'v', heights = c(1,2))
  
  return(PlotFinal)
  
}

MakeGSEA_HeatmapSubsetGO <- function(tf, TargetTerms){
  
  df <- fgsea_GRN_CEN[[tf]]
  df <- subset(df, padj<=0.1)
  colnames(df)[1] <- "Network"
  
  GOs_DB_tem <- unique(subset(GOs_DB, .id == tf)[,c(2,3,6,8)])
  df <- left_join(df, unique(GOs_DB_tem[,c(1,3)]), by=c("pathway"="GO.GRN"))  
  df <- left_join(df, unique(GOs_DB_tem[,c(2,4)]), by=c("pathway"="GO.CEN"))
  
  df$Sig.Terms.GRN[is.na(df$Sig.Terms.GRN)] <- 0
  df$Sig.Terms.CEN[is.na(df$Sig.Terms.CEN)] <- 0
  
  
  df <- subset(df, Sig.Terms.GRN > 0 & Sig.Terms.CEN > 0)
  
  Dic <- unique(df[, c("pathway","Term")])
  
  
  
  print(paste0("... ", tf, " ..."))
  
  if( length(unique(df$Network)) >1 ){
    
    df <- df %>%
      select(Network, NES, pathway) 
    # dcast to calculate order
    # M_DF <- data.table::dcast(df[,1:3], Term ~ Network, value.var = "NES")
    M_DF <- as_tibble(as.data.frame(matrix(0, nrow = length(unique(df$pathway)), ncol = 1)))
    M_DF[,1] <- unique(df$pathway)
    colnames(M_DF)[1] <- "pathway"
    c=2
    nets <- unique(df$Network)
    print(length(nets))
    for (n in nets){
      #
      tem <- subset(df, Network==n)[,-c(1)]
      M_DF <- left_join(M_DF, tem, by="pathway")
      colnames(M_DF)[c] <- n
      print(paste0("... n: ", c-1, " ..."))
      print(dim(M_DF))
      c=c+1
    }
    
    M_DF[is.na(M_DF)] <- 0
    print(dim(M_DF))
    df <- gather(M_DF, value, key, -pathway)
    colnames(df) <- c("pathway", "Network", "NES")
    
    print(dim(df))
    df <- left_join(df, Dic, by="pathway")
    print(dim(df))
    
    
  }
  
  else{
    df <- df %>%
      select(Network, NES, Term) 
  }
  
  
  
  print(paste0(" . Ranking ", tf))
  df <- df %>%
    group_by(Network) %>%
    mutate(RankGO=rank(-NES)) %>%
    arrange(RankGO)  
  
  Rank_df <- df %>% 
    select(Term, RankGO) %>%
    group_by(Term) %>%
    summarise(RankT=median(RankGO)) %>%
    arrange(RankT)
  
  Top <- Rank_df$Term#[1:50]
  #Tail <- Rank_df$Term[(nrow(Rank_df)-20):nrow(Rank_df)]
  #print(Top),  #print(Tail)
  
  DFhot <- subset(df, Term %in% c(Top))
  #DFhot$Term <- factor(DFhot$Term, levels = c(Top, Tail))
  
  net_with_gsea <- unique(DFhot$Network)
  Exp <- TFs_Exp[[tf]]
  Exp <- subset(Exp, Net %in% net_with_gsea) # Print only expression from data expression
  
  Order_Exp <- Exp %>%
    group_by(Net) %>%
    mutate(Mean=mean(CPM)) %>%
    select(Net, Mean) %>%
    arrange(desc(Mean))  
  
  Order_Exp <- unique(Order_Exp)
  Order_Exp <- subset(Order_Exp, Mean >=0.1)$Net
  
  Exp$Net <- factor(Exp$Net, levels = Order_Exp)
  
  Plot1<- ggplot(Exp, aes(x=Net, y=CPM))+
    geom_jitter(size=0.2, alpha=0.3, width = 0.01)+
    theme_pubclean() +
    #scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), position = "right") +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1, size = 8),
          plot.subtitle = element_text(size=8)) +
    xlab("Network") +
    labs(subtitle = paste0(tf, " (", ReplaceName(tf), ")"))
  
  
  DFhot$NES[DFhot$NES > 3] <- 3
  DFhot$NES[DFhot$NES < -3] <- -3
  
  DFhot$Term <- factor(DFhot$Term, levels = c(rev(Rank_df$Term)))
  DFhot <- subset(DFhot, Network %in% Order_Exp) # heatmap if pass exp filter 
  DFhot$Network <- factor(DFhot$Network, levels = Order_Exp)
  print(DFhot)
  
  #DFhot <- filter(DFhot, grepl(paste(TargetTerms, collapse="|"), pathway)) 
  
  Plot2 <- ggplot(DFhot, aes(x=Network, y=Term, fill=NES)) + 
    geom_tile() +
    geom_text(aes(label=RankGO), size=1) +
    scale_fill_viridis_c(direction = 1, option = "A", limits=c(-3,3)) +
    theme_bw()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 7)) 
  
  PlotFinal <- Plot1+ Plot2 +  plot_layout(heights = c(1, 5))
  #PlotFinal <- ggarrange(Plot1, Plot2, 
  #                      nrow = 2, align = 'v', heights = c(1,2))
  
  return(PlotFinal)
  
}

MakeGSEA_HeatmapSubsetInverse <- function(tf, TargetTerms){
  
  df <- fgsea_GRN_CEN[[tf]]
  df <- subset(df, padj<=0.1)
  colnames(df)[1] <- "Network"
  
  GOs_DB_tem <- unique(subset(GOs_DB, .id == tf)[,c(2,3,6,8)])  
  df <- left_join(df, unique(GOs_DB_tem[,c(1,3)]), by=c("pathway"="GO.GRN"))  
  df <- left_join(df, unique(GOs_DB_tem[,c(2,4)]), by=c("pathway"="GO.CEN"))
  
  df$Sig.Terms.GRN[is.na(df$Sig.Terms.GRN)] <- 0
  df$Sig.Terms.CEN[is.na(df$Sig.Terms.CEN)] <- 0
  
  
  df <- subset(df, Sig.Terms.GRN > 0 & Sig.Terms.CEN > 0)
  
  Dic <- unique(df[, c("pathway","Term")])
  print(paste0("... ", tf, " ..."))
  
  
  if( length(unique(df$Network)) >1 ){
    
    df <- df %>%
      select(Network, NES, pathway) 
    # dcast to calculate order
    # M_DF <- data.table::dcast(df[,1:3], Term ~ Network, value.var = "NES")
    M_DF <- as_tibble(as.data.frame(matrix(0, nrow = length(unique(df$pathway)), ncol = 1)))
    M_DF[,1] <- unique(df$pathway)
    colnames(M_DF)[1] <- "pathway"
    c=2
    nets <- unique(df$Network)
    print(length(nets))
    for (n in nets){
      #
      tem <- subset(df, Network==n)[,-c(1)]
      M_DF <- left_join(M_DF, tem, by="pathway")
      colnames(M_DF)[c] <- n
      print(paste0("... n: ", c-1, " ..."))
      print(dim(M_DF))
      c=c+1
    }
    
    M_DF[is.na(M_DF)] <- 0
    print(dim(M_DF))
    df <- gather(M_DF, value, key, -pathway)
    colnames(df) <- c("pathway", "Network", "NES")
    
    print(dim(df))
    df <- left_join(df, Dic, by="pathway")
    print(dim(df))
    
    
  }
  
  else{
    df <- df %>%
      select(Network, NES, Term) 
  }
  
  print(paste0(" . Ranking ", tf))
  df <- df %>%
    group_by(Network) %>%
    mutate(RankGO=rank(-NES)) %>%
    arrange(RankGO)  
  
  Rank_df <- df %>% 
    select(Term, RankGO) %>%
    group_by(Term) %>%
    summarise(RankT=median(RankGO)) %>%
    arrange(RankT)
  
  Top <- Rank_df$Term 
  DFhot <- subset(df, Term %in% c(Top))
  
  net_with_gsea <- unique(DFhot$Network)
  Exp <- TFs_Exp[[tf]]
  Exp <- subset(Exp, Net %in% net_with_gsea) # Print only expression from data expression
  
  Order_Exp <- Exp %>%
    group_by(Net) %>%
    mutate(Mean=mean(CPM)) %>%
    select(Net, Mean) %>%
    arrange(desc(Mean))  
  
  Order_Exp <- unique(Order_Exp)
  Order_Exp <- subset(Order_Exp, Mean >=0.1)$Net
  
  Exp$Net <- factor(Exp$Net, levels = Order_Exp)
  
  Plot1<- ggplot(Exp, aes(x=Net, y=CPM))+
    geom_jitter(size=0.2, alpha=0.3, width = 0.01)+
    theme_pubclean() +
    #scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0), position = "right") +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1, size = 8),
          plot.subtitle = element_text(size=8)) +
    xlab("Network") +
    labs(subtitle = paste0(tf, " (", ReplaceName(tf), ")"))
  
  
  DFhot$NES[DFhot$NES > 3] <- 3
  DFhot$NES[DFhot$NES < -3] <- -3
  
  DFhot$Term <- factor(DFhot$Term, levels = c(rev(Rank_df$Term)))
  
  # Heatmap if pass exp filter 
  DFhot <- subset(DFhot, Network %in% Order_Exp)                
  DFhot$Network <- factor(DFhot$Network, levels = Order_Exp)
  print(DFhot)
  
  DFhot <- filter(DFhot, !grepl(paste(TargetTerms, collapse="|"), Term)) ### 
  
  Plot2 <- ggplot(DFhot, aes(x=Network, y=Term, fill=NES)) + 
    geom_tile() +
    geom_text(aes(label=RankGO), size=1) +
    scale_fill_viridis_c(direction = 1, option = "A", limits=c(-3,3)) +
    theme_bw()+
    theme(panel.border = element_blank(), 
          panel.grid.major = element_blank(), 
          axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 7)) 
  
  PlotFinal <- Plot1+ Plot2 +  plot_layout(heights = c(1, 5))
  #PlotFinal <- ggarrange(Plot1, Plot2, 
  #                      nrow = 2, align = 'v', heights = c(1,2))
  
  return(PlotFinal)
  
}

# ARF7 (Zm00001d039267): phenylpropanoid 
# Based on visual inspection 


gsea_MetabolicTop <- as_tibble(read.table("Candidate_Metabolic_Function_by_TF_GSEA.txt", h=T, sep = "\t"))

TopTerms <- left_join(gsea_MetabolicTop, TFdic, by=c("TF"="TF.Name"))
TopTerms <- unique(TopTerms[,c("TF.v4", "Term")])
TopTerms <- split(TopTerms$Term, TopTerms$TF.v4, )

ReplaceName(names(fgsea_GRN_CEN))

MakeGSEA_HeatmapSubset("Zm00001d031717", TopTerms[["Zm00001d031717"]])

for (tf in names(TopTerms)){
  MakeGSEA_HeatmapSubset(tf, TopTerms[[tf]])
  ggsave(paste0("Figures_/Met_GSEA_",ReplaceName(tf),"_",tf,"",".pdf"), width=8, height=6, units = c("in"))
}


## Tables examples with only top and tail 2 PWY by TF and CoExp Net

MakeGSEA_TopTail <- function(tf){
  
  df <- fgsea_GRN_CEN[[tf]]
  df <- subset(df, padj<=0.1)
  colnames(df)[1] <- "Network"
  
  dfsize <- unique(df[,c("Network", "pathway", "size")])
  
  # do not print GOs which does not have target in both networks
  GOs_DB_tem <- unique(subset(GOs_DB, .id == tf)[,c(2,3,6,8)])
  
  df <- left_join(df, unique(GOs_DB_tem[,c(1,3)]), by=c("pathway"="GO.GRN"))  
  df <- left_join(df, unique(GOs_DB_tem[,c(2,4)]), by=c("pathway"="GO.CEN"))
  
  df$Sig.Terms.GRN[is.na(df$Sig.Terms.GRN)] <- 0
  df$Sig.Terms.CEN[is.na(df$Sig.Terms.CEN)] <- 0
  
  df <- subset(df, Sig.Terms.GRN > 0 & Sig.Terms.CEN > 0)
  
  Dic <- unique(df[, c("pathway","Term")])
  
  print(paste0("... ", tf, " ..."))
  
  print(".... 1 ....")
  print(df)
  if( length(unique(df$Network)) >1 ){
    
    df <- df %>%
      select(Network, NES, pathway) 
    # dcast to calculate order
    # M_DF <- data.table::dcast(df[,1:3], Term ~ Network, value.var = "NES")
    M_DF <- as_tibble(as.data.frame(matrix(0, nrow = length(unique(df$pathway)), ncol = 1)))
    M_DF[,1] <- unique(df$pathway)
    colnames(M_DF)[1] <- "pathway"
    c=2
    nets <- unique(df$Network)
    print(length(nets))
    for (n in nets){
      #
      tem <- subset(df, Network==n)[,-c(1)]
      M_DF <- left_join(M_DF, tem, by="pathway")
      colnames(M_DF)[c] <- n
      print(paste0("... n: ", c-1, " ..."))
      print(dim(M_DF))
      c=c+1
    }
    
    M_DF[is.na(M_DF)] <- 0
    
    
    print(".... 2 ....")
    # print(M_DF)
    df <- gather(M_DF, value, key, -pathway)
    colnames(df) <- c("pathway", "Network", "NES")
    
    print(dim(df))
    df <- left_join(df, Dic, by="pathway")
    print(dim(df))
    print(".... 3 ....")
    # print(df)
    
    
  }
  
  else{
    df <- df %>%
      select(Network, NES, Term) 
  }
  
  print(paste0(" . Ranking ", tf))
  
  df <- subset(df, NES != 0)
  
  df <- df %>%
    group_by(Network) %>%
    mutate(RankGO=rank(-NES, ties.method = 'first')) %>%
    arrange(RankGO)
  
  print(".... 4 ....")
  #print(df)

  print(".... 5 ....")
  
  df_top <- df %>%
    group_by(Network) %>%
    arrange(RankGO) %>%
    slice_head(n=2)
    
  df_tail <- df %>%
    group_by(Network) %>%
    arrange(RankGO) %>%
    slice_tail(n=2)
  
  df_hot <- rbind(df_top, df_tail)
  
  df_hot <- left_join(df_hot, unique(GOs_DB_tem[,c(1,3)]), by=c("pathway"="GO.GRN"))  
  df_hot <- left_join(df_hot, unique(GOs_DB_tem[,c(2,4)]), by=c("pathway"="GO.CEN"))
  colnames(df_hot)[c(6,7)] <- c("Targ_GRN", "Targ_CEN")
  
  df_hot <- left_join(df_hot, dfsize, by=c("Network", "pathway")) 
  return(df_hot)
  
  
}

MYB31_TopTail <- MakeGSEA_TopTail("Zm00001d006236") # MYB31
View(MYB31_TopTail)

##############################################################
#################      Manual pathway ########################
##############################################################

TopTerms <- list()
# ARF7
TopTerms[["Zm00001d039267"]] <-  c("phenylpropanoid", "wall", "flavonoid", "glucuronoxylan", "xylan", "coumarin", "lignin")
# bHLH47
TopTerms[["Zm00001d034298"]] <-  c("monosaccharide", "hexose", "serine", "aromatic", "disaccharide", "isoprenoid", "fatty", "pyruvate", "starch", "wall","disaccharide", "mannose")
# bHLH91
TopTerms[["Zm00001d047017"]] <-  c("chitin", "jasmonic", "salicylic", "ethylene", "abscisic")
# bHLH145 Zm00001d031717
TopTerms[["Zm00001d031717"]] <-  c("jasmonic", "ethylene")
# GSEA_bZIP113_Zm00001d026398
TopTerms[["Zm00001d026398"]] <-  c("phenylpropanoid", "flavonoid", "lignin")
# GSEA_C2H2133 Zm00001d016793
TopTerms[["Zm00001d016793"]] <-  c("jasmonic", "ethylene")
# COL2 Zm00001d033719:  photosynthesis, carotenoid
TopTerms[["Zm00001d033719"]] <-  c("photosynthesis", "carotenoid", "tetraterpenoid", "glyceraldehyde", "isopentenyl")
# GSEA_COL7 Zm00001d025770
TopTerms[["Zm00001d025770"]] <-  c("photosynthesis", "carotenoid", "tetraterpenoid", "glucose", "glyceraldehyde", "isopentenyl", "serine")
# GSEA_COL8 Zm00001d013443
TopTerms[["Zm00001d013443"]] <-  c("photosynthesis", "carotenoid", "tetraterpenoid", "glucose", "glyceraldehyde", "isopentenyl", "serine")
# GSEA COL13 Zm00001d046925
TopTerms[["Zm00001d046925"]] <-  c("flavonoid", "phenylpropanoid", "photosynthesis")
# GSEA_G2 Zm00001d039260: starch, FA,  glucosinolate and glycosinolate 
TopTerms[["Zm00001d039260"]] <-  c("fatty", "glycosinolate", "glucosinolate", "glycosinolate", "starch")
# GSEA_GLK1 Zm00001d044785: penthouse and glucosinolate related processes, terpenoid metabolism
TopTerms[["Zm00001d044785"]] <-  c("penthouse", "glucosinolate", "terpenoid", "glycosinolate", "pentose", "pentose", "pyruvate", "glucose")
# GSEA_GLK52 Zm00001d026542: phenylpropanoid, cysteine, serine, ion homeostasis
TopTerms[["Zm00001d026542"]] <-  c("phenylpropanoid", "cysteine", "serine", "homeostasis", "pigment", "serine")
# GSEA_HB70 Zm00001d025964
TopTerms[["Zm00001d025964"]] <-  c("glyceraldehyde", "NAPD", "glucose", "pentose", "shunt", "Anthocyanin", "phenylpropanoid")
# GSEA_KNOX6 Zm00001d015549
TopTerms[["Zm00001d015549"]] <-  c("phenylpropanoid", "lipid", "flavonoid", "isoprenoid")
# GSEA_MYB6 Zm00001d041576: phenylpropanoid, carotenoid, and FA
TopTerms[["Zm00001d041576"]] <-  c("phenylpropanoid", "carotenoid", "flavonoid", "fatty", "lipid")


MakeGSEA_HeatmapSubset("Zm00001d015549", TopTerms[["Zm00001d015549"]])
MakeGSEA_HeatmapSubsetInverse("Zm00001d015468", TopTerms[["Zm00001d015468"]])
RankedTable <- MakeTableRankingGOs(TF2test[43])


MakeGSEA_HeatmapSubset("Zm00001d039267", TopTerms[["Zm00001d039267"]])
MakeGSEA_HeatmapSubset("Zm00001d034298", TopTerms[["Zm00001d034298"]])
MakeGSEA_HeatmapSubset("Zm00001d047017", TopTerms[["Zm00001d047017"]])
MakeGSEA_HeatmapSubset("Zm00001d031717", TopTerms[["Zm00001d031717"]])
MakeGSEA_HeatmapSubset("Zm00001d026398", TopTerms[["Zm00001d026398"]])
MakeGSEA_HeatmapSubset("Zm00001d016793", TopTerms[["Zm00001d016793"]])
MakeGSEA_HeatmapSubset("Zm00001d033719", TopTerms[["Zm00001d033719"]])
MakeGSEA_HeatmapSubset("Zm00001d025770", TopTerms[["Zm00001d025770"]])

MakeGSEA_HeatmapSubset("Zm00001d046925", TopTerms[["Zm00001d046925"]])
MakeGSEA_HeatmapSubset("Zm00001d039260", TopTerms[["Zm00001d039260"]])
MakeGSEA_HeatmapSubset("Zm00001d044785", TopTerms[["Zm00001d044785"]])
MakeGSEA_HeatmapSubset("Zm00001d026542", TopTerms[["Zm00001d026542"]])
MakeGSEA_HeatmapSubset("Zm00001d025964", TopTerms[["Zm00001d025964"]])
MakeGSEA_HeatmapSubset("Zm00001d015549", TopTerms[["Zm00001d015549"]])
MakeGSEA_HeatmapSubset("Zm00001d041576", TopTerms[["Zm00001d041576"]])

TF2test_Top <- c("Zm00001d039267", "Zm00001d034298", "Zm00001d047017",
                 "Zm00001d031717", "Zm00001d026398", "Zm00001d016793",
                 "Zm00001d033719", "Zm00001d025770", "Zm00001d013443",
                 "Zm00001d046925","Zm00001d039260", "Zm00001d044785",
                 "Zm00001d026542", "Zm00001d025964", "Zm00001d015549",
                 "Zm00001d041576")
##############################################################



