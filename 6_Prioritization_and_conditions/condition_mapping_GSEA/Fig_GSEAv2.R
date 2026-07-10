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
library(pheatmap)
library(topGO)
library(ggpointdensity)
library(factoextra)
library(FactoMineR)
library(dendextend)
library(ggdark)
library(corrplot)

set.seed(42)

###################################################################################
#######                        Functions                                    #######
###################################################################################

ReplaceName <- function(ids){
  
  for (i in 1:nrow(TFdic)){
    ids <- gsub(TFdic$V2[i], TFdic$V1[i], ids)
  }
  return(ids)
}

ReplaceNetIndex <- function(ids){
  
  ids <- CoExpNetAnno$Type[CoExpNetAnno$nid == ids]
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

## 
chop=function(myStr,mySep,myField){
  ## chop a string by a separator and return specified field
  choppedString=sapply(strsplit(myStr,mySep),"[",myField)
  if(length(myField)>1){
    choppedString=apply(choppedString,2,function(x){paste0(x[!is.na(x)],collapse=mySep)})
  }
  return(choppedString)
}


# Calculaye GFA by TFs using wPCC DB from CoExpTFdb   
Get.GSEA.byTF <- function(tf){
  
  # list all TFs wPCC fileson DB
  CoexDB <- list.files(path = "CoExpTFdb/", pattern = tf)
  
  # keep names of files for tf
  CoexDB <- CoexDB[grepl(tf, CoexDB)]
  
  # Get GOs associated with TF 
  GO_targ <-  subset(TFGOs_net, TF == tf)
  
  # GSEA for each GO terms
  GenesInGO <- subset(GeneGO_NetDB, pGO.ID %in% GO_targ$GO)
  
  # List of target genes by GO
  target.list <- split(GenesInGO$GeneID, GenesInGO$pGO.ID)
  
  ## Defined Cor files to read
  file_pcc <-  paste("CoExpTFdb/", CoexDB, sep = "")
  
  # Read co-expression files
  file_pcc <- lapply(file_pcc, function(x) 
                                fread(x, h=T, stringsAsFactors = F)[,2:3]
                                  %>% dplyr::arrange(-wPCC))
  
  # Add names to list of wPCC files 
  names(file_pcc) <- lapply(CoexDB, function(x) chop(x, '[.]', 2))
  
  # Set rank input
  Rank <- lapply(file_pcc, function(x) setNames(x$wPCC, x$GeneID))
  print("... Done reading .. ")
  
  # Fgsea for list of PCC ranks 
  fgseaRes_n <- lapply(Rank, function(x) 
                              fgseaMultilevel(pathways = target.list,
                              stats= x,
                              minSize  = 5, 
                              maxSize  = 2000,
                              eps=0)[,1:7])
  
  # Combine all results by TF
  fgseaRes_n <- rbindlist(fgseaRes_n, idcol = T)
  
  # order by P values
  fgseaRes_n <- fgseaRes_n %>% 
    dplyr::group_by(.id) %>%
    dplyr::arrange(pval, .by_group = TRUE)
  
  fgseaRes_n[,'TF'] <- tf
  write.table(fgseaRes_n, paste0("GSEA_results/GSEA_GOs.", tf, ".txt"), quote = F, sep = "\t", row.names = F)
  
  print(paste0("... Done GSEA: ",tf," .."))

}

Plot.GSEA.byTF <- function(tf, CoExpID, GOid){
  #CoExpID = "n17a_2"
  # list all TFs wPCC fileson DB
  #  wPCC.n16b.Zm00001d042202.txt
  CoexDB <- paste0("CoExpTFdb/wPCC.", CoExpID, ".", tf, ".txt")
  print(" .. step 1 ..")
  # GSEA for each GO terms
  GenesInGO <- subset(GeneGO_NetDB, pGO.ID %in% GOid)
  print(" .. step 2 ..")
  # List of target genes by GO
  target.list <- split(GenesInGO$GeneID, GenesInGO$pGO.ID)
  print(" .. step 3 ..")
  
  
  # Read co-expression files
  file_pcc <- fread(CoexDB, h=T, stringsAsFactors = F)[,2:3] %>% dplyr::arrange(-wPCC)
  print(" .. step 4 ..")
  # Set rank input
  Rank <- setNames(file_pcc$wPCC, file_pcc$GeneID)
  print("... Done reading .. ")
  
  # Fgsea for list of PCC ranks 
  fgseaRes_n <- fgseaMultilevel(pathways = target.list,
                                stats= Rank,
                                minSize  = 5, 
                                maxSize  = 2000,
                                eps=0)[,1:7]
  
  
  Plot <- plotEnrichment(target.list[[1]], Rank) + labs(title=subset(TFGOs_net, TF == tf & GO == GOid)$term)
  
  print(paste0("... Done GSEA: ",tf," .."))
  return(Plot)
}

MakeGSEA_Heatmap  <- function(tf) {
  # Used GSEAdb
  #tf=Phenyl_related
  temGSEAdb <- subset(GSEAdb, TF == tf)
  
  # reduce the NES value in not significant association
  temGSEAdb$NES[is.na(temGSEAdb$NES)]  <- 0
  rm(tf)
  #temGSEAdb$NES[temGSEAdb$padj > 0.1] <- temGSEAdb$NES[temGSEAdb$padj > 0.1]*0.1
  
  #
  temGSEAdb <- left_join(temGSEAdb, GOdic, by= c('pathway'='parent'))[,c(".id",'NES','parentTerm','pathway')]
  temGSEAdb <- temGSEAdb[!(is.na(temGSEAdb$parentTerm)),]
  
  # Matrix without scaling to draw heatmap with original values
  M.temGSEAdb <- dcast(temGSEAdb, parentTerm ~ .id, value.var = 'NES') 
  row.names(M.temGSEAdb) <- M.temGSEAdb$parentTerm
  M.temGSEAdb <- M.temGSEAdb[, -c(1)]
  M.temGSEAdb[is.na(M.temGSEAdb)] <- 0
  cat('.. Done Matrix ..')
  
  colnames(M.temGSEAdb) <- as.character(sapply(colnames(M.temGSEAdb), ReplaceNetIndex)) 
  cat('.. Done net index ..')
  
  # # Add TFs names as column
  # MyAnn <- data.frame(Net.Type=CoExpNetAnno[colnames(M.temGSEAdb),'Type'])
  # row.names(MyAnn) <- colnames(M.temGSEAdb)
  
  
  Plot_h <- pheatmap(M.temGSEAdb,
                     #annotation_row=CoExpNetAnno[,4], 
                     #annotation_col = MyAnn,
                     cluster_rows = T,
                     #cutree_rows =3,
                     fontsize_col = 5,
                     cellheight = 8,
                     cellwidth= 6, 
                     show_rownames=T,
                     show_colnames=T,
                     #annotation_colors = my_colour,
                     color = viridis(n = 100, option = 'D', direction = 1, alpha = 0.8),
                     treeheight_col=0,
                     treeheight_row=0,
                     border_color = NA,
                     scale = "none")
  return(Plot_h)
  
}

MakeGSEA_CorPlot  <- function(tf) {
  # Used GSEAdb
  #tf=ABA_related
  temGSEAdb <- subset(GSEAdb, TF == tf)
  
  # reduce the NES value in not significant association
  temGSEAdb$NES[is.na(temGSEAdb$NES)]  <- 0
  
  #temGSEAdb$NES[temGSEAdb$padj > 0.1] <- temGSEAdb$NES[temGSEAdb$padj > 0.1]*0.1
  
  #
  temGSEAdb <- left_join(temGSEAdb, GOdic, by= c('pathway'='parent'))[,c(".id",'NES','parentTerm','pathway')]
  temGSEAdb <- temGSEAdb[!(is.na(temGSEAdb$parentTerm)),]
  
  # Remaple parenTerm
  
  temGSEAdb[,"parentTerm2"] <- paste0("[",temGSEAdb$pathway, "]", " ", temGSEAdb$parentTerm)
  
  # Matrix without scaling to draw heatmap with original values
  M.temGSEAdb <- dcast(temGSEAdb, parentTerm2 ~ .id, value.var = 'NES') 
  row.names(M.temGSEAdb) <- M.temGSEAdb$parentTerm2
  M.temGSEAdb <- M.temGSEAdb[, -c(1)]
  M.temGSEAdb[is.na(M.temGSEAdb)] <- 0
  cat('.. Done Matrix ..')
  
  colnames(M.temGSEAdb) <- as.character(sapply(colnames(M.temGSEAdb), ReplaceNetIndex)) 
  cat('.. Done net index ..')
  
  
  # # Add TFs names as column
  # MyAnn <- data.frame(Net.Type=CoExpNetAnno[colnames(M.temGSEAdb),'Type'])
  # row.names(MyAnn) <- colnames(M.temGSEAdb)
  
  CorGOs <- cor(t(M.temGSEAdb))
  CorNets <- cor(M.temGSEAdb)
  
  
  col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                             "cyan", "#007FFF", "blue", "#00007F"))
  
  Plot_GOs <- {corrplot(CorGOs,  order = "hclust", addrect = 2, col = col1(100), tl.cex=0.7)}
  Plot_Nets <- corrplot(CorNets,  order = "hclust", addrect = 4, col = col1(100), tl.cex=0.5)
  
  name_1 <- paste0('Plots/Fig_FigS15.', ReplaceName(tf),'.GOs.pdf')
  name_2 <- paste0('Plots/Fig_FigS15.', ReplaceName(tf),'.Nets.pdf')
  
  pdf(file = name_1)
  corrplot(CorGOs, type = "lower", mar = c(0,0,1,0),  number.cex = 0.5, tl.col = "black",
           number.digits = 2, order = "hclust",  addrect = 2, col = col1(100),
           tl.cex=0.7, pch.col = 'black') 
  dev.off()
  
  pdf(file = name_2)
  corrplot(CorNets, type = "lower", mar = c(0,0,1,0),  number.cex = 0.5, tl.col = "black",
           number.digits = 2, order = "hclust",  addrect = 2, col = col1(100),
           tl.cex=0.7, pch.col = 'black') 
  dev.off()
  
  
  
}

Plot.kmeanSize.byTF <- function(tf, ExpID, GOid){
  # ExpID = "n17a_2"
  # GOid = "GO:0006635"
  # tf = Lipid_related
  # list all TFs wPCC fileson DB
  # wPCC.n16b.Zm00001d042202.txt
  expDBid <- paste0("ExpDB/", ExpID, ".tsv")
  print(" .. step 1 ..")
  
  # GSEA for each GO terms
  GenesInGO <- unique(subset(GeneGO_NetDB, pGO.ID %in% GOid))
  print(" .. step 2 ..")
  
  # List of target genes by GO
  target.list <- unique(GenesInGO$GeneID)
  print(" .. step 3 ..")
  
  
  # Read co-expression files
  expDB <- read.table(expDBid, h=T, stringsAsFactors = F)
  gids <- expDB$gid
  expDB <- expDB[,-c(1)]
  row.names(expDB) <- gids
  
  # Expression filter
  mask <- (rowSums((expDB >= 0.5)*1) >= ncol(expDB)*0.1)
  expDB <- expDB[mask,]
  gids <- names(mask)[mask == TRUE]
  
  # scale expression
  expDB <- scale(expDB)
  target.list <- target.list[target.list %in% row.names(expDB)]
  expDB <- expDB[c(tf, target.list),]
  #dim(expDB)
  #expDB <- subset(expDB, gid %in% c(tf))
  print(" .. step 4 ..")
  
  # Explore possible number of cluster: Visualize eigenvalues/variances
  nclust <- fviz_nbclust(expDB, kmeans, method = "gap_stat")
  return(nclust)
}

Top.SamplesExpression.byTF <- function(tf, ExpID){
  # ExpID = "n14e"
  # tf = Lipid_related
  
  # list all TFs wPCC fileson DB
  # wPCC.n16b.Zm00001d042202.txt
  expDBid <- paste0("ExpDB/", ExpID, ".tsv")
  print(" .. step 1 ..")
  
  # Read co-expression files
  expDB <- read.table(expDBid, h=T, stringsAsFactors = F)
  gids <- expDB$gid
  expDB <- expDB[,-c(1)]
  row.names(expDB) <- gids
  
  
  # Expression filter
  mask <- (rowSums((expDB >= 0.5)*1) >= ncol(expDB)*0.1)
  expDB <- expDB[mask,]
  gids <- names(mask)[mask == TRUE]
  
  # slace  expression
  expDB <- expDB <- scale(expDB)
  expDB <- expDB[c(tf),]
  print(" .. step 2 ..")
  
  # reshape DF
  longDF <- reshape2::melt(as.matrix(expDB)) %>% 
    as_tibble() %>%
    dplyr::rename('GeneID'='Var2') %>%
    dplyr::rename('Var2'="Var1")
  
  longDF$GeneID <- tf
  
  ExpAnnDB <- as.data.table(read.table("ExpressionSamples_annotation.txt", sep = "\t", h=T))
  
  # Add Expression ann to gene df
  longDF$Var2 <- as.character(longDF$Var2)
  
  longDF <- left_join(longDF, ExpAnnDB, by=c('Var2'="SampleID"))
  
  
  longDF$GeneID <- as.character(longDF$GeneID)
  longDF$Var2 <- as.character(longDF$Var2)
  
  # defined a single sample name
  longDF$Cluster <- paste0(longDF$Tissue, ":", longDF$Genotype, ":", longDF$Treatment)
  
  # sort samples by TF Expression
  subset(longDF, GeneID==tf) %>%
    dplyr::arrange(value) %>%
    dplyr::select(Var2) -> SampleClass
  
  SampleClass[,'SampleID'] <- seq(1, nrow(SampleClass), 1)
  
  #
  longDF <- left_join(longDF, SampleClass, by=c("Var2"))
  longDF$Var2 <- factor(longDF$Var2, levels = SampleClass$Var2)
  
  longDF %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(GeneID=unique(GeneID),
                     Var2=unique(Var2),
                     value=mean(value),
                     SampleID=min(SampleID)) -> longDF
  
  longDF %>%
    dplyr::arrange(value) %>%
    dplyr::top_n(5, -value)  -> Taildf
  
  longDF %>%
    dplyr::arrange(value) %>%
    dplyr::top_n(5, value)  -> Topdf
  
  longDF <- rbind(Topdf, Taildf)
  
  return(longDF)
}

Top.SamplesExpression.byTF(Lipid_related, "n19a")

Plot.Expression.byTF <- function(tf, ExpID, GOid, k, Add_Label){
  # ExpID = "n19a"
  # GOid = "GO:0006635"
  # tf = Lipid_related
  # k=4
  
  # list all TFs wPCC fileson DB
  # wPCC.n16b.Zm00001d042202.txt
  expDBid <- paste0("ExpDB/", ExpID, ".tsv")
  print(" .. step 1 ..")
  
  
  # GSEA for each GO terms
  GenesInGO <- unique(subset(GeneGO_NetDB, pGO.ID %in% GOid))
  print(" .. step 2 ..")
  
  
  # List of target genes by GO
  target.list <- unique(GenesInGO$GeneID)
  print(" .. step 3 ..")
  
  # Read co-expression files
  expDB <- read.table(expDBid, h=T, stringsAsFactors = F)
  gids <- expDB$gid
  expDB <- expDB[,-c(1)]
  row.names(expDB) <- gids
  
  # Expression filter
  mask <- (rowSums((expDB >= 0.5)*1) >= ncol(expDB)*0.1)
  expDB <- expDB[mask,]
  gids <- names(mask)[mask == TRUE]
  
  # scale expression
  expDB <- scale(expDB)
  
  # defined genes in go which are also expressed in dataset
  target.list <- target.list[target.list %in% row.names(expDB)]
  
  #Subset genes of interes
  expDB <- expDB[c(tf, target.list),]
  
  #expDB <- subset(expDB, gid %in% c(tf))
  print(" .. step 4 ..")
  
  # Explore possible number of cluster: Visualize eigenvalues/variances
  km.res <- kmeans(expDB, k, nstart = 25)
  km.res <- as.data.frame(km.res$cluster)
  km.res[,'GeneID'] <- row.names(km.res)
  
  #
  colnames(km.res)[1] <- c("Cluster")
  print(as.data.table(table(km.res$Cluster)))
  
  longDF <- reshape2::melt(as.matrix(expDB)) %>% 
    as_tibble() %>%
    dplyr::rename('GeneID'="Var1") %>%
    dplyr::rename('Var2'="Var2")
  
  # add cluster to gene df
  longDF <- left_join(longDF, km.res, by='GeneID')
  longDF$GeneID <- as.character(longDF$GeneID)
  longDF$Var2 <- as.character(longDF$Var2)
  
  #
  longDF$Cluster <- paste0('C', longDF$Cluster)
  
  # sort samples by TF Expression
  subset(longDF, GeneID==tf) %>%
    dplyr::arrange(value) %>%
    dplyr::select(Var2) -> SampleClass
  
  SampleClass[,'SampleID'] <- seq(1, nrow(SampleClass), 1)
  
  #
  longDF <- left_join(longDF, SampleClass, by=c("Var2"))
  longDF$Var2 <- factor(longDF$Var2, levels = SampleClass$Var2)
  
  longDF %>%
    dplyr::select(Cluster, value) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::summarise(Mean=mean(value), Sd=sd(value)) -> DFsummary
    
  print(DFsummary)
    
  
  TopDF <- Top.SamplesExpression.byTF(tf, ExpID) 
  
  #
  tittlename <- paste0(subset(TFGOs_net, TF == tf & GO == GOid)$term,"; ",
                       subset(CoExpNetAnno, nid == ExpID)$Type)
  
  #
  if (Add_Label==T) {
    ggplot(longDF, aes(x=SampleID, y=value, color=Cluster)) +
      geom_line(data=subset(longDF, GeneID==tf), 
                aes(x=SampleID, y=value), color='black', linetype = "dashed") +
      geom_smooth(fill = "grey", linewidth=0.5)  +
      scale_color_brewer(palette = "Set1") +
      geom_point(data = TopDF, aes(x=SampleID, y=value), color = "red") +
      geom_text_repel(data = TopDF, aes(x=SampleID, y=value, label=Cluster),
                      color='Black', size = 2,
                      box.padding = 0.5, max.overlaps = Inf) +
      ylab('Scale CPM') +  xlab("Samples") +
      #theme_pubclean() +
      theme_pubclean()+
      #dark_theme_minimal() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) + 
      labs(title=tittlename) +
      theme(strip.text.x = element_text(size = 10), 
            axis.text=element_text(size=10),
            plot.subtitle=element_text(size=3, hjust=0, color="black"),
            text = element_text(size=10),
            legend.position="bottom",
            legend.key.size = unit(0.3, 'cm')) -> Plotout
    
  }
  
  else {
    ggplot(longDF, aes(x=SampleID, y=value, color=Cluster)) +
      geom_line(data=subset(longDF, GeneID==tf), 
                aes(x=SampleID, y=value), color='black', linetype = "dashed") +
      geom_smooth(fill = "grey", linewidth=0.5)  +
      scale_color_brewer(palette = "Set1") +
      ylab('Scale CPM') +  xlab("Samples") +
      #theme_pubclean() +
      theme_pubclean()+
      #dark_theme_minimal() +
      scale_x_continuous(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0)) + 
      labs(title=tittlename) +
      theme(strip.text.x = element_text(size = 10), 
            axis.text=element_text(size=10),
            plot.subtitle=element_text(size=3, hjust=0, color="black"),
            text = element_text(size=10),
            legend.position="bottom",
            legend.key.size = unit(0.3, 'cm')) -> Plotout
  }
  
  return(Plotout)
}


###################################################################################

##############################################################################
##################         Read data input           #########################
##############################################################################

# TF Names
TFdic <- as_tibble(read.table("Data/Annotations/TF_Id_Name.txt", h=F, stringsAsFactors = F))
colnames(TFdic)

# GOs term annotations, already reduced/ mapped to parents as used on Fig_pecanpyPart2
GeneGO_NetDB <- as_tibble(readRDS('../Fig_PecanpyPart2/MaizeSyntenicGenes_GOparent.rds'))

# TFs-GO network to test
TFGOs_net <- as_tibble(readRDS('../Fig_PecanpyPart2/TFGO_4337_net.rds')) 
TF2test <- unique(TFGOs_net$TF)

# TF names
GOdic <- unique(fread("../Fig_MethodsComparison/ReduceGOterms_All_methods.txt",  header =T)[,1:2])

write.table(data.frame(GeneID=TF2test),
            '/maindisk/fabio/Projects/MaizeENCODE/Data_45_net/wPCC_net_only_TFs/TF_list.txt',
            row.names = F, sep = '\t', quote = F)

CoExpNetAnno <- read.table('CoExpNet_annotation.txt', sep = '\t', header = T)
row.names(CoExpNetAnno) <- CoExpNetAnno$net_index
CoExpNetAnno[,'Type'] <- paste0(CoExpNetAnno$net_index, "::", CoExpNetAnno$net_type, "::", CoExpNetAnno$note)

##############################################################################

################################################################
########      Calculate GSEA score for TF-GO net         ####### 
################################################################

# Nu
tested <- list.files(path = "GSEA_results/", pattern = '^GSEA_*')
tested <- gsub('GSEA_GOs.', '', tested)
tested <- unique(gsub('.txt', '', tested))

# Check if TFs already tested
GeneTF2test <- TF2test
GeneTF2test <- TF2test[!(TF2test %in% tested)]
#

Lgenes <- length(GeneTF2test) # Maximum rank: value used in loop
w = 3# Size of range to test

print(".. Ready to start ..")
for (i in seq(0, Lgenes, w)){
  max=Lgenes
  Start=i+1
  end=i+w
  if (end<max){
    listtotest <- GeneTF2test[Start:end]
    print(length(listtotest))
    Range <- paste(Start,end, sep = "-")
    # list of genes to test and number of cores
    mclapply(listtotest, Get.GSEA.byTF, mc.cores=w)
    print(paste("... Done gsea ", Range, " ...."))
  }
  else if (end>max) {
    break
  } 
  else{
    r <- paste(Start, max, sep = "-")
    listtotest <- GeneTF2test[Start:max]
    print(length(listtotest))
    Range <- paste(Start,max, sep = "-")
    w = max - Start
    # list of genes to test and number of cores
    mclapply(listtotest, wCorRows, mc.cores=w)
    print(paste("... Done gsea ", Range, " ...."))
  }
}

################################################################

################################################################
########         Count GOs enriched using GSEA          ######## 
################################################################

GSEAdb <- lapply(paste0("GSEA_results/", list.files(path = "GSEA_results/", pattern = '^GSEA_*')), 
                 fread)

# Cobine all results 
GSEAdb <- rbindlist(GSEAdb)

##################################
#### Totals Nets tested by TF ####
##################################

NetsByTF <- GSEAdb %>% 
  dplyr::select(.id, TF) %>%
  unique() %>%
  dplyr::select(TF) %>%
  table() %>%
  as.data.table() %>%
  dplyr::rename(TotalNets=N) %>%
  dplyr::rename(TF=".") 

# Totals Nets tested by TF with at least a significan hit
NetsByTF_sign <- GSEAdb %>% 
  dplyr::filter(padj <= 0.1) %>%
  dplyr::select(.id, TF) %>%
  unique() %>%
  dplyr::select(TF) %>%
  table() %>%
  as.data.table() %>%
  dplyr::rename(TotalNetsSig=N) %>%
  dplyr::rename(TF=".") 

NetsByTF <- left_join(NetsByTF, NetsByTF_sign, by='TF') %>%
                dplyr::mutate(Ratio=round((TotalNetsSig/TotalNets)*100, 2))

NetsByTF
##################################

##################################
#### Totals GO tested by TF   ####
##################################

GOsByTF <- GSEAdb %>% 
  dplyr::select(.id, pathway, TF) %>% # Selected GO, net, and TF
  unique() %>%                        # keep unique values
  dplyr::select(.id, TF) %>%          # Use net (.id) and TF to counts Freq of GOs (total)
  table() %>%          
  as.data.table() %>%
  dplyr::filter(N > 0) %>%            # Keep net and TF tested at least once (N > 0)
  dplyr::rename(TotalGOs=N) %>%
  dplyr::rename(Net=.id) 

GOsByTF_sign <- GSEAdb %>% 
  dplyr::filter(padj <= 0.1) %>%
  dplyr::select(.id, pathway, TF) %>% # Selected GO, net, and TF
  unique() %>%                        # keep unique values
  dplyr::select(.id, TF) %>%          # Use net (.id) and TF to counts Freq of GOs (total)
  table() %>%          
  as.data.table() %>%
  dplyr::filter(N > 0) %>%            # Keep net and TF tested at least once (N > 0)
  dplyr::rename(SigGOs=N) %>%
  dplyr::rename(Net=.id) 

##
GOsByTF <- left_join(GOsByTF, GOsByTF_sign, by=c('Net', 'TF')) %>%
  dplyr::mutate(Ratio=round((SigGOs/TotalGOs)*100, 2))

GOsByTF$Ratio[is.na(GOsByTF$Ratio)] <- 0

GOsByTF %>%
  dplyr::select(TF, Ratio) %>%
  dplyr::group_by(TF) %>%
  dplyr::summarise(GOmean=mean(Ratio)) -> Mean_GOsByTF

mean(Mean_GOsByTF$GOmean)



##################################

#################################################
####  Cluster TFs bases on their GSEA results ###
#################################################
# https://rpkgs.datanovia.com/factoextra/

# make wide matrix with % of GOs terms by TF and Net
M.GOsByTF <- dcast(GOsByTF[,c(1,2,5)], TF ~ Net, value.var = 'Ratio') 
row.names(M.GOsByTF) <- M.GOsByTF$TF
M.GOsByTF <- M.GOsByTF[, -c(1)]
M.GOsByTF[is.na(M.GOsByTF)] <- 0

# Scale values to better identification of cluster of TF
M.GOsByTF <- scale(M.GOsByTF)

# Explore possible number of cluster: Visualize eigenvalues/variances
fviz_nbclust(M.GOsByTF, kmeans, method = "gap_stat")
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
res.pca <- prcomp(iris[, -5],  scale = TRUE)


km.res <- kmeans(M.GOsByTF, 3, nstart = 25)

fviz_cluster(km.res, data = M.GOsByTF,
             palette = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07", "#9370DB"), 
             labelsize = 3,
             ggtheme = theme_minimal(),
             main = "Partitioning Clustering Plot")

# Visualize
fviz_dend(res, rect = TRUE, cex = 0.2, 
          color_labels_by_k = FALSE,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800", "#FC4E07", "9370DB"))


# Compute hierarchical clustering and cut into 4 clusters
res <- hcut(M.GOsByTF, k = 3, stand = TRUE)
Cluster <- cutree(tree = as.dendrogram(res), k=3)
Cluster <- as.data.frame(Cluster)
Cluster$Cluster <- paste0('Cluster', Cluster$Cluster)

# Rename cluster to mach order in heatmap
Cluster$Cluster <- gsub('Cluster2', 'Cluster4', Cluster$Cluster)
Cluster$Cluster <- gsub('Cluster3', 'Cluster2', Cluster$Cluster)
Cluster$Cluster <- gsub('Cluster4', 'Cluster3', Cluster$Cluster)
table(Cluster$Cluster)

# Matrix without scaling to draw heatmap with original values
M2.GOsByTF <- dcast(GOsByTF[,c(1,2,5)], TF ~ Net, value.var = 'Ratio') 
row.names(M2.GOsByTF) <- M2.GOsByTF$TF
M2.GOsByTF <- M2.GOsByTF[, -c(1)]
M2.GOsByTF[is.na(M2.GOsByTF)] <- 0

# Add TFs names as column
Cluster[,'TF'] = rownames(Cluster) 
Cluster <- as.data.table(Cluster)
#################################################

#######################################################
####  Map number of GO by TF into clusters to        ##
####  test bias by number of GOs and cluster classes ##
#######################################################

# Count GOs by TF
GOsByTF_red <- GSEAdb %>% 
  dplyr::select(pathway, TF) %>%  # Selected GO, net, and TF
  unique() %>%                         # keep unique values
  dplyr::select(TF) %>%           # Use net (.id) and TF to counts Freq of GOs (total)
  table() %>%          
  as.data.table() %>%
  dplyr::filter(N > 0) %>%             # Keep net and TF tested at least once (N > 0)
  dplyr::rename(TotalGOs=N) %>%
  dplyr::rename(TF=".") 

# Count GOs by TF
Cluster <- left_join(Cluster, GOsByTF_red, by='TF')

#######################################################

########################################################
####  Test differences in average percentage of GO   ###
####  significant by TF                              ###
#### (n=total nets testes. Val=percentage by net).   ###
########################################################

left_join(Mean_GOsByTF, Cluster, by="TF") %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(Mean=mean(GOmean))

left_join(left_join(Mean_GOsByTF, Cluster, by="TF"), GOsByTF[,c("Net", "TF", "Ratio")], by="TF") %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(Mean=mean(Ratio))

#######################################################


######################################################## 
####                     Figures                    ####
########################################################

########
#### Plot 1: Percentage of Nets with at least a Sig GO
########
NetsByTF$Ratio[is.na(NetsByTF$Ratio)] <- 0

NetsByTF %>%
  ggplot(aes(x="TFs", y=Ratio)) +
  geom_boxplot(notch = T, outlier.shape = NA, fill='#6666FF', alpha=0.5) +
  #geom_jitter(alpha=0.7, size=0.1, width = 0.1) +
  theme_pubclean() +
  #ylab(bquote(GO[Z~";TFs per GO"])) + 
  ylab("Percentage Coexp. nets\nwith >= 1 Sig. GO term") + 
  xlab("") +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        plot.subtitle=element_text(size=6, hjust=0, color="black"),
        text = element_text(size=10)) -> Plot_1

mean(NetsByTF$Ratio)
########

########
#### Plot 2: Average percentage of GO significant by TF (n=total nets testes. Val=percentage by net)
########
Mean_GOsByTF %>%
  ggplot(aes(x="TFs", y=GOmean)) +
  geom_boxplot(notch = T, outlier.shape = NA, fill='#E5FFCC', alpha=0.5) +
  #geom_jitter(alpha=0.7, size=0.1, width = 0.1) +
  theme_pubclean() +
  #ylab(bquote(GO[Z~";TFs per GO"])) + 
  ylab("Average Per. of Sig.\nGO terms by TF") + 
  xlab("") +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        plot.subtitle=element_text(size=6, hjust=0, color="black"),
        text = element_text(size=10)) -> Plot_2

##
Mean_GOsByTF %>%
  ggplot(aes(x=TotalGOs, y=GOmean)) +
  #geom_point()
  geom_pointdensity() +
  scale_color_viridis() +
  theme_pubclean() +
  #ylab(bquote(GO[Z~";TFs per GO"])) +
  ylab("") +
  xlab("Total GO terms ber TF") +
  theme(strip.text.x = element_text(size = 10),
        axis.text=element_text(size=10),
        plot.subtitle=element_text(size=6, hjust=0, color="black"),
        text = element_text(size=10)) -> Plot_2b


# stats
Mean_GOsByTF %>%
  dplyr::filter(GOmean <= quantile(Mean_GOsByTF$GOmean, 0.25) & GOmean > 0 ) %>%
  dplyr::arrange(GOmean)

Mean_GOsByTF %>%
  dplyr::filter(GOmean >= quantile(Mean_GOsByTF$GOmean, 0.75) ) %>%
  dplyr::arrange(GOmean)

#GOsByTF[GOsByTF$TF == 'Zm00001d021389',]
#TFGOs_net[TFGOs_net$TF == 'Zm00001d021389',]

########

########
#### Plot 3: Heatmap with cluster of genes
########
my_colour = list(Cluster = c(Cluster1 = "#FF69B4", Cluster2 = "#00CED1", Cluster3='#CC99FF'))

Plot_3 <- pheatmap(M2.GOsByTF,
         annotation_row=Cluster[,1], 
         cluster_rows = res,
         cutree_rows =3,
         cellwidth=2.5, 
         cellheight = 0.3,
         show_rownames=F,
         show_colnames=F,
         annotation_colors = my_colour,
         color = viridis(n = 100, option = 'B', direction = 1, alpha = 0.8),
         treeheight_col=0,
         treeheight_row=0, 
         scale = "none")
########

########
#### Plot 4: GOs by TF between clusters
########

my_comparisons <- list( c("Cluster1", "Cluster2"), c("Cluster2", "Cluster3"), c("Cluster1", "Cluster3"))

ggplot(Cluster, aes(x=Cluster, y=TotalGOs, fill=Cluster)) +
  geom_boxplot(notch = T, outlier.size = 0.5) +
  #geom_jitter(alpha=0.7, size=0.1, width = 0.1) +
  theme_pubclean() +
  scale_fill_manual(values = c(Cluster1 = "#FF69B4", Cluster2 = "#CC99FF", Cluster3='#00CED1')) +
  #ylab(bquote(GO[Z~";TFs per GO"])) + 
  ylab("GOs by TF") + 
  xlab("") +
  stat_compare_means(comparisons=my_comparisons, 
                     label="p.signif", 
                     method='wilcox.test') +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        plot.subtitle=element_text(size=6, hjust=0, color="black"),
        text = element_text(size=10), 
        legend.position = "none") -> Plot_4

Plot_4
########

########
#### Plot 5: Sig. GOs by TF
########
ClusterDF <- Cluster
ClusterDF[,"TF"] <- rownames(ClusterDF)

# Significant GOs by TF by cluster
left_join(Mean_GOsByTF, ClusterDF, by="TF") %>%
  ggplot(aes(x=Cluster, y=GOmean, fill=Cluster)) +
  geom_boxplot(notch = T, outlier.size = 0.5) +
  scale_fill_manual(values = c(Cluster1 = "#FF69B4", Cluster2 = "#CC99FF", Cluster3='#00CED1')) +
  theme_pubclean() +
  #ylab(bquote(GO[Z~";TFs per GO"])) + 
  ylab("Average Per. of Sig.\nGO terms by TF") + 
  xlab("") +
  stat_compare_means(comparisons=my_comparisons,  label="p.signif", method='wilcox.test') +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        plot.subtitle=element_text(size=6, hjust=0, color="black"),
        text = element_text(size=10), 
        legend.position = "none") -> Plot_5

########

########
#### Plot 6: Sig. GOs by TF between clusters
########

# Add cluster into to TFs and NES information
left_join(NetsByTF, ClusterDF, by="TF") %>%
  ggplot(aes(x=Cluster, y=TotalNetsSig, fill=Cluster)) +
  geom_boxplot(notch = T, outlier.size = 0.5) +
  scale_fill_manual(values = c(Cluster1 = "#FF69B4", Cluster2 = "#CC99FF", Cluster3='#00CED1')) +
  theme_pubclean() +
  # ylab(bquote(GO[Z~";TFs per GO"])) + 
  ylab("Percentage Network with\nSig. GO terms by TF") + 
  xlab("") +
  stat_compare_means(comparisons=my_comparisons,  label="p.signif", method='wilcox.test') +
  theme(strip.text.x = element_text(size = 10), 
        axis.text=element_text(size=10),
        plot.subtitle=element_text(size=6, hjust=0, color="black"),
        text = element_text(size=10), 
        legend.position = "none") -> Plot_6

Plot_6

left_join(NetsByTF, ClusterDF, by="TF") %>%
  dplyr::group_by(Cluster) %>%
  dplyr::filter(is.na(TotalNetsSig)==FALSE) %>%
  dplyr::summarise(Avg=mean(TotalNetsSig))
########

########
#### Plot 7: heatmap with TF examples 
########


## Defined TF of interest for further analysis
subset(TFGOs_net, FDR <= 0.1) 

ABA_related    <- TFdic[TFdic$V1 %in% c('bHLH43'),]$V2
Lipid_related  <- TFdic[TFdic$V1 %in% c('ARF14'),]$V2 # PRH115, 'Zm00001d052815'
Phenyl_related <- TFdic[TFdic$V1 %in% c('HB33'),]$V2
Leaf_related   <- TFdic[TFdic$V1 %in% c('WRKY25'),]$V2

targetTFs <- TFdic[TFdic$V1 %in% c("bHLH43", "ARF14", "HB33", "WRKY25"),]$V2


# bHLH43 Zm00001d033267
# ARF14  Zm00001d050781
# WRKY25 Zm00001d032265
# HB33   Zm00001d033378

GOsByTF[GOsByTF$TF == 'Zm00001d032265',]
TFGOs_net[TFGOs_net$TF == 'Zm00001d050781',]

subset(GSEAdb, TF == 'Zm00001d033267' & padj <= 0.05)

ClusterDF[targetTFs, ]

# Subset targeted TFs
M2_targetTFs <- M2.GOsByTF[targetTFs,]
colnames(M2_targetTFs) <- as.character(sapply(colnames(M2_targetTFs), ReplaceNetIndex)) 


Plot_7 <- pheatmap(M2_targetTFs,
                   annotation_row=ClusterDF[ClusterDF$TF %in% targetTFs,], 
                   cluster_rows = TRUE,
                   cutree_rows =3,
                   cellwidth=8, 
                   cellheight = 10,
                   show_rownames=F,
                   show_colnames=T,
                   annotation_colors = my_colour,
                   color = viridis(n = 100, option = 'B', direction = 1, alpha = 0.8),
                   treeheight_col=0,
                   treeheight_row=0, 
                   scale = "none")
Plot_7


########
#### Plot 8: GSEA specific examples
########

dev.off()
Plot_ABA <- MakeGSEA_Heatmap(ABA_related)       # type 1
Plot_Lipid <- MakeGSEA_Heatmap(Lipid_related)   # type 1
Plot_Phenyl <- MakeGSEA_Heatmap(Phenyl_related) # type 1
Plot_Leaf <- MakeGSEA_Heatmap(Leaf_related)     # type 2

# Corplot based on NES 
MakeGSEA_CorPlot(ABA_related)         # type 1
MakeGSEA_CorPlot(Lipid_related)       # type 1
MakeGSEA_CorPlot(Phenyl_related)       # type 1
MakeGSEA_CorPlot(Leaf_related)       # type 1

tem <- subset(GSEAdb, TF == Lipid_related)
tem <- left_join(tem, GOdic, by= c('pathway'='parent'))[,c(".id",'NES','parentTerm','pathway')]
tem[,"Index"] <- as.character(sapply(tem$.id, ReplaceNetIndex)) 
View(tem)


########


########
#### Plot 9: Expression profiles for examples of ARF14  specific examples
########

left_join(NetsByTF, ClusterDF, by="TF") %>%
  dplyr::arrange(-TotalNetsSig) 


TFGOs_net[TFGOs_net$TF %in% Lipid_related,]

# GO:0006635 lipids
# GO:0009699 Phenolics
# GO:0006569 tryptophan 

left_join(subset(GSEAdb, TF==Lipid_related) %>% 
            filter(pathway=="GO:0006635") %>%
            filter(NES >= max(NES)*0.5 | NES <= -1.3), 
          CoExpNetAnno[,c('Type', "net_index", "nid")], by =c(".id"="nid")) %>%
  dplyr::arrange(-NES)

left_join(subset(GSEAdb, TF==Lipid_related) %>% 
            filter(pathway=="GO:0006569") %>%
            filter(NES==max(NES) | NES <= -1.3), 
          CoExpNetAnno[,c('Type', "net_index", "nid")], by =c(".id"="nid")) %>%
  dplyr::arrange(-NES)

left_join(subset(GSEAdb, TF==Lipid_related) %>% 
            filter(pathway=="GO:0006569") %>%
            filter(NES==max(NES) | NES <= -1.3), 
          CoExpNetAnno[,c('Type', "net_index", "nid")], by =c(".id"="nid")) %>%
  dplyr::arrange(-NES)


# Plot: ARF14
Plot_Lipid_related_Top <- Plot.GSEA.byTF(Lipid_related, 'n18e_2', 'GO:0006635') # N28
Plot_Phe_related_Top <- Plot.GSEA.byTF(Lipid_related, 'n18e_2', 'GO:0009699')   # N28
Plot_Tryp_related_Top <- Plot.GSEA.byTF(Lipid_related, 'n18e_2', 'GO:0006569')  # N28

Plot_Lipid_related_tail <- Plot.GSEA.byTF(Lipid_related, 'n18e_1', 'GO:0006635') # N27
Plot_Phe_related_tail <- Plot.GSEA.byTF(Lipid_related, 'n18e_1', 'GO:0009699')   # N27
Plot_Tryp_related_tail <- Plot.GSEA.byTF(Lipid_related, 'n18e_1', 'GO:0006569')  # N27

Plot_Lipid_related_tail2 <- Plot.GSEA.byTF(Lipid_related, 'n19a', 'GO:0006635') # N5 N11
Plot_Phe_related_tail2 <- Plot.GSEA.byTF(Lipid_related, 'n19a', 'GO:0009699')   # N5
Plot_Tryp_related_tail2 <- Plot.GSEA.byTF(Lipid_related, 'n19a', 'GO:0006569')  # N5

Top.SamplesExpression.byTF(Lipid_related, "n19a")

Plot_9a  <- Plot_Lipid_related_Top/Plot_Phe_related_Top/Plot_Tryp_related_Top
Plot_9b  <- Plot_Lipid_related_tail/Plot_Phe_related_tail/Plot_Tryp_related_tail
Plot_9c  <- Plot_Lipid_related_tail2/Plot_Phe_related_tail2/Plot_Tryp_related_tail2

Plot_9  <-  { Plot_9a | Plot_9b | Plot_9c }

########
#### Plot 10: Expression profiles for examples of ARF14  specific examples
########
#
Plot.kmeanSize.byTF(Lipid_related, 'n18e_2', 'GO:0006635') # N28:Lipids
Plot.kmeanSize.byTF(Lipid_related, 'n18e_2', 'GO:0009699') # N28
Plot.kmeanSize.byTF(Lipid_related, 'n18e_2', 'GO:0006569') # N28
#
Plot_Lipid_T1_Exp_AFR14 <-  Plot.Expression.byTF(Lipid_related, 'n18e_2', 'GO:0006635', 7, TRUE)
Plot_Phe_T1_Exp_AFR14   <-  Plot.Expression.byTF(Lipid_related, 'n18e_2', 'GO:0006635', 4, FALSE)
Plot_Tryp_T1_Exp_AFR14  <-  Plot.Expression.byTF(Lipid_related, 'n18e_2', 'GO:0006635', 4, FALSE)

Plot.kmeanSize.byTF(Lipid_related, 'n18e_1', 'GO:0006635') # N27:Lipids
Plot.kmeanSize.byTF(Lipid_related, 'n18e_1', 'GO:0009699') # N27
Plot.kmeanSize.byTF(Lipid_related, 'n18e_1', 'GO:0006569') # N27
#
Plot_Lipid_T1_Exp_AFR14 <-  Plot.Expression.byTF(Lipid_related, 'n18e_1', 'GO:0006635', 3, TRUE)
Plot_Phe_T1_Exp_AFR14   <-  Plot.Expression.byTF(Lipid_related, 'n18e_2', 'GO:0006635', 4, FALSE)
Plot_Tryp_T1_Exp_AFR14  <-  Plot.Expression.byTF(Lipid_related, 'n18e_2', 'GO:0006635', 4, FALSE)

#
Plot.kmeanSize.byTF(Lipid_related, 'n19a', 'GO:0006635') # N3:Lipids
Plot.kmeanSize.byTF(Lipid_related, 'n19a', 'GO:0009699') # N3:Phe
Plot.kmeanSize.byTF(Lipid_related, 'n19a', 'GO:0006569') # N3:Tryp

Plot_Lipid_T2_Exp_AFR14 <-  Plot.Expression.byTF(Lipid_related, 'n19a', 'GO:0006635', 4, TRUE)
Plot_Phe_T2_Exp_AFR14   <-  Plot.Expression.byTF(Lipid_related, 'n19a', 'GO:0009699', 3, FALSE)
Plot_Tryp_T2_Exp_AFR14  <-  Plot.Expression.byTF(Lipid_related, 'n19a', 'GO:0006569', 3, FALSE)

Plot10 <- (Plot_Lipid_T2_Exp_AFR14 / Plot_Phe_T2_Exp_AFR14 / Plot_Tryp_T2_Exp_AFR14) & theme(legend.position = 'right') # + plot_layout(guides = "collect")
Plot10

########


########################################################


##################################

##############
# Save plots
##############


Fig_FigS13 <- (Plot_1 | Plot_2 | Plot_2b) + plot_layout(widths = c(1, 1, 4), guides = "collect") & theme(legend.position = 'bottom')
Fig_FigS13def <-  Plot_5 | Plot_6 | Plot_4


pdf("Plots/Fig_FigS13.pdf", width=7, height=4)
print(Fig_FigS13)
dev.off()

pdf("Plots/Fig_FigS13def.pdf", width=8, height=3)
print(Fig_FigS13def)
dev.off()

pdf("Plots/Fig_Fig4a.pdf", width=8, height=5)
print(Plot_3)
dev.off()

pdf("Plots/Fig_Fig4b.pdf", width=8, height=5)
print(Plot_7)
dev.off()

pdf("Plots/Fig_Fig4c.pdf", width=7, height=4)
print(Plot_Lipid)
dev.off()

pdf("Plots/Fig_Fig4d.pdf", width=4, height=5)
print(Plot_9)
dev.off()

pdf("Plots/Fig_Fig4e.pdf", width=4, height=5)
print(Plot10)
dev.off()

#
pdf("Plots/Fig_FigS14a.pdf", width=7, height=5)
print(Plot_ABA)
dev.off()

pdf("Plots/Fig_FigS14b.pdf", width=7, height=5)
print(Plot_Phenyl)
dev.off()

pdf("Plots/Fig_FigS14c.pdf", width=6, height=4)
print(Plot_Leaf)
dev.off()

#############


##################################

# Total GOs that pass in at least a what! 
TotalGOs <- as.data.table(table(TFGOs_net$TF)) %>%
  dplyr::rename(Total=N) %>%
  dplyr::rename(TF=V1) 

# Count number of times a GO was tested by TF
GO_timesTested <- as.data.table(table(GSEAdb[,c(2,9)])) %>% 
  dplyr::filter(N>0) %>%
  dplyr::rename(Total=N)


# Calculate freq of GO sig. enrichment in CoExNet by TF
GONetgsea_Freq <- GSEAdb[,c(1,2,4,9)]

# padj as binary
GONetgsea_Freq$padj <- (GONetgsea_Freq$padj <= 0.1)*1

# Count freq
GONetgsea_Freq <- as.data.table(table(GONetgsea_Freq[,2:4])) %>% 
                    dplyr::filter(N>0)

# add total times GO tested
GONetgsea_Freq <- left_join(GONetgsea_Freq, GO_timesTested, 
                            by = c('pathway', 'TF'))

# 
GOsByTF




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

MakeGSEA_HeatmapSubset("Zm00001d015549", TopTerms[["Zm00001d015549"]])
MakeGSEA_HeatmapSubsetInverse("Zm00001d015468", TopTerms[["Zm00001d015468"]])

MakeGSEA_HeatmapSubset("Zm00001d039267", TopTerms[["Zm00001d039267"]])

TF2test_Top <- c("Zm00001d039267", "Zm00001d034298", "Zm00001d047017",
                 "Zm00001d031717", "Zm00001d026398", "Zm00001d016793",
                 "Zm00001d033719", "Zm00001d025770", "Zm00001d013443",
                 "Zm00001d046925","Zm00001d039260", "Zm00001d044785",
                 "Zm00001d026542", "Zm00001d025964", "Zm00001d015549",
                 "Zm00001d041576")

##############################################################

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
  
  df <- df %>%
    group_by(Network) %>%
    mutate(RankGO=rank(-NES)) %>%
    arrange(RankGO)
  
  
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





######
# Size: 20x15
MakeGSEA_Heatmap(TF2test[6])
MakeGSEA_Heatmap(TF2test[38])
MakeGSEA_Heatmap(TF2test[43])
